#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /cvmfs/icecube.opensciencegrid.org/users/gen2-optical-sim/software/icetray/build

from os.path import expandvars
from icecube import icetray, dataclasses, dataio, phys_services, trigger_sim
from icecube.icetray import OMKey
from icecube.dataclasses import I3OMGeo, I3ModuleGeo
from icecube import DOMLauncher
from icecube import vuvuzela   
from icecube import sim_services
from icecube.gen2_sim.utils import ModuleToString, BasicHitFilter, iceprod_wrapper
from icecube.gen2_sim.segments.trigger import RunUpgradeTriggers
try:
    from icecube import mDOM_WOM_simulation
except: 
    print("Unable to import mDOM_WOM_simulation, probably because clsim client/server stuff screwed up the calls")
    pass


from I3Tray import I3Units
import numpy as np
# from icecube.icetray import I3Units, logging
# logging.set_level(
#     logging.I3LogLevel.LOG_DEBUG)

'''
The input pulse map to the function mcpulse_to_recopulse need to be sorted in time in order for it to bin pulses correctly. 
This following function sorts I3MCPulseMap which is then called as an input for the mcpulse_to_recopulse function
 -- NLad 
'''

def sort_pulsemap(unsorted_map):
    sorted_map = type(unsorted_map)()
    for omkey, pulses in unsorted_map:
        sorted_map[omkey] = type(pulses)(sorted(pulses, key=lambda pulse:pulse.time))
    return sorted_map

def mcpulse_to_recopulse(frame, mapname = "I3MCPulseSeriesMapGen2", outputmap = "I3RecoPulseSeriesMapGen2",bin_width=1*I3Units.nanosecond,):
    '''
        A module that does a direct conversion of I3MCPulses to I3RecoPulses.
        It is intended to use on the PMTResponseSimulation output when you
        when one wants to avoid the DOM simulation for some reason (no DOM electronic simulation. ie no launches but
        PMT effects such as saturation is present).
    '''
    from icecube.icetray import I3Units
    recopulsemap = dataclasses.I3RecoPulseSeriesMap()
    mcpulsemap = frame[mapname]
    sorted_mcpulsemap = sort_pulsemap(mcpulsemap)
    for omkey, pulses in sorted_mcpulsemap:
        # In order for some more advanced reconstructions to work
        # properly the pulses need to be coaleced
        if(len(pulses)==0):
            continue
        recopulsemap[omkey] = dataclasses.I3RecoPulseSeries()
        c = pulses[0].charge
        t = pulses[0].time*c
        last_t = pulses[0].time
        for p in pulses[1:]:
            
            #Adding pulse to bin if it is within the bin_width time
            if( (p.time - last_t) < bin_width):
                c += p.charge
                t += p.time*p.charge
                last_t = t/c
            else:
                #creating a new recopulse of the coaleced pulses
                rpulse = dataclasses.I3RecoPulse()
                rpulse.time = t/c
                rpulse.charge = c
                        
                rpulse.width = bin_width
                rpulse.flags = dataclasses.I3RecoPulse.PulseFlags.LC
                if rpulse.charge > 0.25:
                    recopulsemap[omkey].append(rpulse)
                c = p.charge
                t = p.time*p.charge
                last_t = t/c
                
        rpulse = dataclasses.I3RecoPulse()
        #same discriminator threshold cut as IC for now
        # if (c>0.25):
        rpulse.time = t/c
        rpulse.charge = c
        rpulse.width = bin_width
        rpulse.flags = dataclasses.I3RecoPulse.PulseFlags.LC
        if rpulse.charge > 0.25:
            recopulsemap[omkey].append(rpulse)
    
    frame[outputmap] = recopulsemap


def scale_noise_rate(frame, scale=1):
    """
    Scale overall noise rate, leaving other parameters untouched
    """
    # NB: If we were to be rigorous and delete I3Calibration before replacing
    # it with a modified copy, the previous version would remain in the
    # previous module's frame-mixing cache for all time. Since these are
    # potentially large (1.3 GB for the mDOM in Sunflower 240!) we cheat and
    # modify in place, secure in the knowledge that no other modules are
    # relying on the constness of the I3Calibration at the point where this is
    # called.
    calib = frame['I3Calibration']
    
    for k in calib.dom_cal.keys():
        if k.string < 87:
            continue
        cal = calib.dom_cal[k]
        cal.dom_noise_thermal_rate *= scale
        cal.dom_noise_decay_rate *= scale
        cal.dom_noise_rate *= scale
        calib.dom_cal[k] = cal


@icetray.traysegment
def AddNoise(tray, name,
             InputName ="I3MCPESeriesMap",
             OutputName = "",
             RandomService = None,
             ExcludeList = [],
             StartTime = -10*I3Units.microsecond,
             EndTime = 10*I3Units.microsecond,
             mDOMGlass = "",
             module = I3OMGeo.mDOM,
             pregeneratedmDOM = False,
             If = lambda f: True):
    """
    Vuvuzela, without Inject (we already have noise parameters)
    """
    # if RandomService == None:
    #     RandomService = phys_services.I3SPRNGRandomService(self.rngseed, self.rngnumberofstreams, self.rngstream)
    # elif type(RandomService) == str:
    RandomService = tray.context[RandomService]

    if module == I3OMGeo.PDOM:
        ThermalRate = 1.73456e-7
        DecayRate = 5.6942e-8
        ScintillationHits = 8.072
        ScintillationMean = 4.395
        ScintillationSigma = 1.777
    elif module == I3OMGeo.mDOM:
        if not mDOMGlass.lower() in ["vitrovex", "benthos"]:
            icetray.logging.log_fatal("Error! You have to specify correct glass type! (vitrovex/benthos)")
        if mDOMGlass.lower()=="vitrovex":
            ThermalRate = 2.07e-8
            DecayRate = 5.0502e-8
            ScintillationHits = 7.55
            ScintillationMean = 2.904
            ScintillationSigma = 1.007
        elif mDOMGlass.lower()=="benthos":   
            ThermalRate = 1.0e-8
            DecayRate = 1.95e-8
            ScintillationHits = 3.91
            ScintillationMean = 3.131
            ScintillationSigma = 1.05
    elif module == I3OMGeo.DEgg:
        ThermalRate = 1.00e-7
        DecayRate = 2.61e-7
        ScintillationHits = 3.6
        ScintillationMean = 7.9
        ScintillationSigma = 3.8
    elif module == I3OMGeo.IceCube:
        ThermalRate = 1.73456e-7
        DecayRate = 5.6942e-8
        ScintillationHits = 8.072
        ScintillationMean = 4.395
        ScintillationSigma = 1.777
        tray.AddModule("Inject", name+"AddNoiseParams",
                       InputNoiseFile = expandvars("$I3_SRC/vuvuzela/resources/data/parameters.dat"),
                   )
        
    # The upgrade modules shouldn't have noise parameters set. Kill whatever's there.
    if module != I3OMGeo.IceCube:
        tray.Add(scale_noise_rate, name+"_killNoiseParams",
                 scale=np.nan, Streams=[icetray.I3Frame.Calibration])

    if (pregeneratedmDOM) and (module == I3OMGeo.mDOM):
        print("WARN: Producing with the experimental new pre-generated mDOM noise.")
        tray.AddModule(vuvuzela.PregeneratedSampler, name+'_mdom_pregenerator',
                       EndWindow = EndTime,
                       InputPath = '/cvmfs/icecube.opensciencegrid.org/users/mlarson/mdom_noise/mDOM_darkrate_-25C.*.npy',
                       ModuleType = I3ModuleGeo.mDOM,
                       OutputHitMap = OutputName, 
                       PhysicsHitMap = InputName,
                       RandomService = RandomService,
                       StartWindow = StartTime, 
                   )
    else:
        tray.AddModule("Vuvuzela", name+"_vuvuzela",
                       InputHitSeriesMapName  = InputName, 
                       OutputHitSeriesMapName = OutputName,
                       StartWindow            = StartTime,
                       EndWindow              = EndTime,
                       RandomService          = RandomService,
                       OMTypes                = [module, ],
                       DOMsToExclude          = ExcludeList,
                       SimulateNewDOMs        = True,
                       DisableLowDTCutoff     = True,
                       UseIndividual          = True,
                       
                       DecayRate = DecayRate,
                       ScintillationHits = ScintillationHits,
                       ScintillationMean = ScintillationMean,
                       ScintillationSigma = ScintillationSigma,
                       ThermalRate = ThermalRate,
                   )



@icetray.traysegment
def DetectorSim(tray, name,
    GCDFile,
    RunID=None,
    RandomServiceName = "I3RandomService",   
    EHEApproximation = False,
    TimeWindow = 0.2,
    KeepMCHits = False,    
    FilterMode= False,
    InputMCPESeriesMaps = ['I3MCPESeriesMap_IceCube',
                           'I3MCPESeriesMap_PDOM',
                           'I3MCPESeriesMap_mDOM',
                           'I3MCPESeriesMap_DEgg',
                       ],
    Sensors = [I3OMGeo.OMType.IceCube,
               I3OMGeo.OMType.PDOM,
               I3OMGeo.OMType.mDOM,
               I3OMGeo.OMType.DEgg,
               ],
    AreaScale = 1.,
    ZeroHitFilter = False,
    VuvuzelaStartTime = -15.*I3Units.microsecond,
    VuvuzelaEndTime = 15.*I3Units.microsecond,
    mDOMGlass="vitrovex",
    pregeneratedmDOM = False,
    ):
    InputMCPESeriesMapName_withNoise = "I3MCPESeriesMapWithNoise"
    upgradeRecoPulses = []
    noisePEList = []
    if ZeroHitFilter:
        # Checking if at least one MCPE was produced, if not - rejecting the event
        # since it will be noise-only 
        tray.Add(BasicHitFilter, 
                 InputMCPESeriesMaps=InputMCPESeriesMaps, 
                 Streams = [icetray.I3Frame.DAQ])
    for i, Sensor in enumerate(Sensors):
        # Adding noise one module type at a time. 
        # No input, since I want to force all of the modules to produce noise in
        # the same window. Doing it the standard way would make the last module type
        # produce noise in +- (interaction window + NModuleTypes * 10 mus), which 
        # would just waste processing time. This explicitly assumes that most events
        # we simulate take less than  10 mus to die off.
        SensorName = ModuleToString(Sensor)
        noise_pe_name = "I3MCPESeriesMapNoise"+SensorName
        noisePEList.extend( [noise_pe_name, InputMCPESeriesMapName_withNoise + SensorName] )

        tray.Add(AddNoise, name+SensorName + "_vuvuzela",
                 InputName = "",
                 OutputName = noise_pe_name,
                 StartTime = VuvuzelaStartTime,
                 EndTime = VuvuzelaEndTime,
                 RandomService = RandomServiceName,
                 mDOMGlass= mDOMGlass,
                 module = Sensor,
                 pregeneratedmDOM = pregeneratedmDOM,
         )
  
        if SensorName not in InputMCPESeriesMaps[i]:
            print(f"WARN: {SensorName} seems mismatched with {InputMCPESeriesMaps[i]}, are they passed in the correct order?")
        # Combine with the physics hits
        tray.Add('I3CombineMCPE', name+SensorName + "_combineMCPE",
                 InputResponses = [InputMCPESeriesMaps[i], noise_pe_name],
                 OutputResponse = InputMCPESeriesMapName_withNoise + SensorName,
                 # If = lambda f: f.Has(InputMCPESeriesMaps[i])
        )

        # This seems to cause a frame collision with InputMCPESeriesMaps[i]
        # If doesn't catch if it exists
        # tray.Add("Copy", name+SensorName+'_copyMCPE',
        #          Keys = [noise_pe_name, InputMCPESeriesMapName_withNoise + SensorName],
        #          If = lambda f: not f.Has(InputMCPESeriesMaps[i]),
        # )
                
        # Merge any MCPEs closer than 200 picoseconds
        tray.AddModule('Delete', 'delete_pidmap'+SensorName,
                       Keys=[InputMCPESeriesMapName_withNoise+SensorName+'ParticleIDMap'])
        tray.AddModule('I3MCPEMerger', name+SensorName+'_coalesceMCPE',
                       Input=InputMCPESeriesMapName_withNoise + SensorName,
                       Output=InputMCPESeriesMapName_withNoise + SensorName,
                       TimeWindow=TimeWindow)
        
        # We now have a series of MCPESeriesMaps. We need to run the
        # detector simulation for each type separately, since its handled 
        # differently for each module

        ####################
        ## Start with the standard IceCube stuff
        ####################
        if Sensor == I3OMGeo.OMType.IceCube:
            tray.Add("PMTResponseSimulator",name+"_PMTResponse_IC86",
                     Input=InputMCPESeriesMapName_withNoise+SensorName,
                     Output="I3MCPulseSeriesMap" + SensorName,
                     MergeHits=True,
                     LowMem = True,
                     EHEApproximation=EHEApproximation,
                     RandomServiceName=RandomServiceName,
                     #If = lambda x: False
            )
            tray.AddModule("DOMLauncher", name+SensorName+"_DOMLauncher",
                           Input= "I3MCPulseSeriesMap" + SensorName,
                           Output="InIceRawData",
                           UseTabulatedPT=True,
                           RandomServiceName=RandomServiceName,
            )

        elif Sensor == I3OMGeo.OMType.PDOM:            
            tray.Add("PMTResponseSimulator",name+"_PMTResponse_PDOM",
                     Input=InputMCPESeriesMapName_withNoise+SensorName,
                     Output="I3MCPulseSeriesMap" + SensorName,
                     MergeHits=True,
                     LowMem = True,
                     EHEApproximation=EHEApproximation,
                     RandomServiceName=RandomServiceName,
                     #SaturationScale = AreaScale,
                     #If = lambda x: False
            )

        elif Sensor == I3OMGeo.OMType.mDOM:
            tray.Add('mDOMPMTSimulator', name+"_PMTSimulator_mDOM",
                     Gain        = 1.0,
                     Input=InputMCPESeriesMapName_withNoise+SensorName,
                     Output="I3MCPulseSeriesMap" + SensorName,
                     TransitTime = 0.0 * I3Units.ns,
                     TTS         =  1.45 * I3Units.ns,
                     ApplySPE    = True
            )
 
        elif Sensor == I3OMGeo.OMType.DEgg:
            # TODO: Is this the right thing to do?
            tray.Add("PMTResponseSimulator","PMTResponse_Gen2degg",
                     Input=InputMCPESeriesMapName_withNoise+SensorName,
                     Output="I3MCPulseSeriesMap" + SensorName,
                     MergeHits = False, 
                     PrePulseProbability = 0.0,
                     LatePulseProbability = 0.0,
                     AfterPulseProbability = 0.0,
                     ApplySaturation = False, 
                     UsePMTJitter = True,
                     EHEApproximation=EHEApproximation,
                     RandomServiceName=RandomServiceName,
                 )
 
        # Only run the simplified module sim for non-IceCube DOMs
        # Can think about using different bin/pulse widths for each module someday
        # by moving the mcpulse_to_recopulse inside of the previous if statements
        if Sensor != I3OMGeo.OMType.IceCube:
            tray.Add(mcpulse_to_recopulse, name+SensorName+"mcpulse_to_recopulse",
                     mapname = "I3MCPulseSeriesMap" + SensorName, 
                     outputmap = "I3RecoPulseSeriesMap" + SensorName,
                     bin_width=10*I3Units.nanosecond, # For now some baseline value
                     Streams = [icetray.I3Frame.DAQ,],
            )
            upgradeRecoPulses.append("I3RecoPulseSeriesMap" + SensorName)

    # Now we have all of the various sensors processed. 
    # Combine the upgrade pulses into one I3RecoPulseSeriesMapGen2 pulsemap
    def merge_pulses(frame, Input=[], Output = 'I3RecoPulseSeriesMapGen2'):
        frame[Output] = dataclasses.I3RecoPulseSeriesMapUnion(frame, Input)
        
    tray.Add(merge_pulses,
             Input = upgradeRecoPulses,
             Streams = [icetray.I3Frame.DAQ])

    # Should we have the triggers in here directly?
    # Yes, temporarily add the triggers here. (BAC 2021-04-14)
    # Though, we should move these into a separate tray segment soon.
    # Be sure to leave OldTriggers = True, so that TriggerSim gets run,
    # and the I3TriggerHierarchy will get created.
    tray.AddSegment(
        RunUpgradeTriggers, name+'_Upgrade_Trigers',
        GCD=dataio.I3File(GCDFile),
        UpgradePulses='I3RecoPulseSeriesMapGen2',
        OldTriggers=True,
        Filter=FilterMode
        )

    # Clean up the photons before we move on
    tray.Add("Delete", Keys = noisePEList + ['I3MCPulseSeriesMap_IceCubeParticleIDMap', 'I3MCPulseSeriesMap_IceCube'])
    if not KeepMCHits:
        tray.AddModule("Delete", name+"_cleanup_I3MCHits_2",
                       Keys = InputMCPESeriesMaps)
    

DetectorSimulation = iceprod_wrapper(DetectorSim)

if __name__ == "__main__":
    stats = {}    
    module = DetectorSimulation()
    module.ExecuteOpts(stats)
    print(stats)
