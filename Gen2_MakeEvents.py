import os,sys
from os.path import expandvars
import logging
import math

from I3Tray import *
from icecube import (icetray, dataio, dataclasses, phys_services,
                     interfaces, simclasses, sim_services)
from icecube.gen2_sim.utils import ModuleToString


class DAQCounter(icetray.I3Module):
    count   = 0
    nevents = 0

    def __init__(self,ctx):
        icetray.I3Module.__init__(self,ctx)
        self.AddParameter("NEvents","name of event counter",self.nevents)
        self.AddOutBox("OutBox");

    def Configure(self):
        self.nevents = self.GetParameter("NEvents")

    def DAQ(self,frame):
        self.count += 1
        if self.count > self.nevents: 
           self.RequestSuspension()
        self.PushFrame(frame)


def BasicCounter(frame,name="basic_counter",Stats={}):
    if name not in Stats:
        Stats[name]=0
    Stats[name] += 1

gcddir = '/cvmfs/icecube.opensciencegrid.org/users/gen2-optical-sim/gcd'
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-g","--gcdfile",type=str,help="GCD-File",required=True,dest="gcdfile")
parser.add_argument("-o","--outputfile",type=str,help="Output file",required=True,dest="outputfile")
# parser.add_argument("-i","--infile",type=str,help="Input File",required=True,dest="infile")

parser.add_argument("--summaryfile",default="",help='XMLSummary filename',dest="summaryfile")
parser.add_argument("--mjd",default="62502",type=str,help='MJD for the GCD file',dest="mjd")
parser.add_argument("--RunId",default=0,type=int,help="Configure run ID",dest="runid")
parser.add_argument("--RNGSeed",default=0,type=int,help="RNG Seed",dest="rngseed")
parser.add_argument("--RNGStream",default=0,type=int,help="RNGStream",dest="rngstream")
parser.add_argument("--RNGNumberOfStreams",default=1,type=int,help="RNGNumberOfStreams",dest="rngnumberofstreams")
parser.add_argument("-n","--nevents",default=50,type=int,help='Number of events',dest="nevents")
parser.add_argument("--SimMode",default="FULL",type=str,help='simulation mode',dest="simmode")
parser.add_argument("--VTXGenMode",default="NuGen",type=str,help='simulation mode',dest="vtxgenmode")
parser.add_argument("--InjectionMode",default="Surface",type=str,help='injection mode',dest="injectionmode")
parser.add_argument("--AutoExtendMuonVolume",default=True,type=bool,help='Detection volume extension, set off for starting events',dest="autoextendmuonvolume")
parser.add_argument("--NuFlavor",required=True,type=str,help='Neutrino Flavor',dest="nuflavor")
parser.add_argument("--gamma",default=2.0,type=float,help='Gamma index',dest="gamma")
parser.add_argument("--FromEnergy",default=1.*I3Units.TeV,type=float,help='Minimum energy',dest="fromenergy")
parser.add_argument("--ToEnergy",default=10.*I3Units.PeV,type=float,help='Maximum energy',dest="toenergy")
parser.add_argument("--ZenithMin",default=0.,type=float,help='min zenith',dest="zenithmin")
parser.add_argument("--ZenithMax",default=180.*I3Units.degree,type=float,help='max zenith',dest="zenithmax")
parser.add_argument("--AzimuthMin",default=0.,type=float,help='min Azimuth',dest="azimuthmin")
parser.add_argument("--AzimuthMax",default=360.*I3Units.degree,type=float,help='max Azimuth',dest="azimuthmax")
parser.add_argument("--CrossSections",default='csms_differential_v1.0',type=str,help='cross section tables',dest="crosssections")
parser.add_argument("--CrossSectionsPath",default=None,help='cross section tables path',dest="crosssectionspath")
parser.add_argument("--ParamsMap",default=dict(),help='any other parameters',dest="paramsmap")
parser.add_argument("--cylinder",default=[0,0,0,0,0],help='For CIRCLE[radius, active_height_before, active_height_after],for SURFACE[radius, length, center_x, center_y, center_z]',dest="cylinderparams")
parser.add_argument("--UseExtrudedInjectionSurface", type=bool,default=True,help='Use MuonGun injection surface extruded from selected GCD file, else default to Cylinder',dest="extrudedsurface")

#parser.add_argument('--gcd', type=str, help='GCDFile to use for simulation',
#                    default=expandvars(
#                        os.path.join(gcddir,
#                                     'IceCubeHEX_Sunflower_240m_v3_generated_IceCubeExtendedDepthRange_mDOM.GCD.i3.bz2')))
parser.add_argument('--gcd1size', type=str, help='GCDFile with equal-size DOMs for photon prop',
                    default=expandvars(
                        os.path.join(gcddir,
                                     'IceCubeHEX_Sunflower_240m_v3_generated_IceCubeExtendedDepthRange.GCD.i3.bz2')))
parser.add_argument('--hybrid', action="store_true", default=False)
parser.add_argument('--sensors', nargs="+", default=['IceCube', 'mDOM'])
parser.add_argument('--effs', type=float, nargs="+", default=[1.,2.2])
parser.add_argument('--no-gpu', dest='gpu', action="store_false", default=True)
parser.add_argument('--dataset', type=int, default=1)
parser.add_argument('--nfiles', type=int, default=1000)
parser.add_argument('--fileno', type=int, default=0)

parser.add_argument('--KeepMCHits', type=bool, default=False,dest="keepmchits")
parser.add_argument('--FilterMode', type=bool, default=False,dest="filtermode")
parser.add_argument('--AreaScale', type=float, default=1.0,dest="areascale")
parser.add_argument('--VuvuzelaStartTime', type=float, default=-15.*I3Units.microsecond,dest="vuvuzelastarttime")
parser.add_argument('--VuvuzelaEndTime', type=float, default=15.*I3Units.microsecond,dest="vuvuzelaendtime")
parser.add_argument('--infile', default=None,dest="infile")
parser.add_argument('--MuonGun', dest='onlyMuonGun', action="store_true", default=False)
parser.add_argument('--UseCRModel',type=bool, default=False, dest="useCRModel")
parser.add_argument("--PowerLawOffset",default=1.0*I3Units.TeV,type=float,help='Power Law Offset for Muon Propogation',dest="ploffset")
parser.add_argument('--SegPath', default="/cvmfs/icecube.opensciencegrid.org/users/gen2-optical-sim/software/icetray/build/lib/icecube/gen2_sim/segments",dest="segpath")

args = parser.parse_args()

sys.path.append(args.segpath)

from icecube.dataclasses import I3OMGeo

from icecube import (neutrino_generator, earthmodel_service, PROPOSAL, cmc, phys_services)
from icecube.simprod.segments import GenerateNeutrinos, PropagateMuons
from icecube.gen2_sim.segments.Calibration import Calibration
from icecube.gen2_sim.segments.BaseProc import BaseProc
from icecube.gen2_sim.segments.BaseReco import BaseReco

from icecube.gen2_sim.segments.clsim import CLSimTraySegment, MakePhotonsMultiSensor, MakePEFromPhotons
# support json ParamsMap, so we can get our dict from the iceprod config file
try:
    import json
except ImportError:
    json = None
if isinstance(args.paramsmap,str):
    if not json:
        raise Exception('ParamsMap provided as string, but python does not understand json')
    args.paramsmap = json.loads(args.paramsmap)

#

# Instantiate a tray 
tray = I3Tray()
icetray.set_log_level(icetray.I3LogLevel.LOG_DEBUG)
randomService = phys_services.I3SPRNGRandomService(args.rngseed, args.rngnumberofstreams, args.rngstream)
tray.context['I3RandomService'] = randomService
tray.context['I3FileStager'] = dataio.get_stagers()


#if(args.infile is None):
#   fnl=[args.gcdfile]
#else:
#   fnl=[args.gcdfile] + [args.infile]


#tray.Add("I3Reader","read", FilenameList = fnl)
# if not args.summaryfile:
#     raise Exception('must define summary file')
# tray.AddService("I3XMLSummaryServiceFactory", "summary",
#     outputfilename=args.summaryfile,
# )
###############muon gun
if(args.onlyMuonGun):
    from icecube.MuonGun import load_model, Floodlight, StaticSurfaceInjector, Cylinder, OffsetPowerLaw, ExtrudedPolygon
    from icecube.MuonGun.segments import GenerateBundles
    from icecube.sim_services import I3ParticleTypePropagatorServiceMap
    from icecube.PROPOSAL import I3PropagatorServicePROPOSAL
    from icecube.cmc import I3CascadeMCService
    # re-use the same RNG for modules that need it on the context
   # tray.context[randomServiceName] = randomService
    def make_propagators(seed):
        propagators = I3ParticleTypePropagatorServiceMap()
        muprop = I3PropagatorServicePROPOSAL()#cylinderHeight=1600.0, cylinderRadius=800.0)
        cprop = I3CascadeMCService(phys_services.I3GSLRandomService(seed)) # dummy RNG
        for pt in 'MuMinus', 'MuPlus':
            propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = muprop
        for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
            propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cprop
        return propagators
    # parameters for MuonGun
    CylinderLength = 1700.0*I3Units.m #FIXME
    CylinderRadius = 1800.0*I3Units.m #FIXME
    CylinderCenter = dataclasses.I3Position(-700.0*I3Units.m, 150.0*I3Units.m, 0.0*I3Units.m)
    PowerLawGamma = args.gamma
    PowerLawOffset = args.ploffset*I3Units.GeV #FIXME
    PowerLawEmin = args.fromenergy*I3Units.GeV #FIXME
    PowerLawEmax = args.toenergy*I3Units.GeV #FIXME




    # Create single muon
    # Create an offset power-law energy distribution of the for dP/dEmu propto (Emu + b)(^-gamma)
    if args.extrudedsurface:
        surface = ExtrudedPolygon.from_file(args.gcdfile, padding=60.0*I3Units.m)
    else:
        surface = Cylinder(CylinderLength, CylinderRadius, CylinderCenter)
    spectrum = OffsetPowerLaw(PowerLawGamma, PowerLawOffset, PowerLawEmin, PowerLawEmax)
    # Set up the generator. This gets stored in a special frame for later reference
    if args.useCRModel:
        model = load_model('GaisserH4a_atmod12_SIBYLL')
        # CR Flux models: Hoerandel5_atmod12_SIBYLL, GaisserH4a_atmod12_SIBYLL, GaisserH4a_atmod12_DPMJET, GaisserH4a_atmod12_DPMJET-C
        model.flux.max_multiplicity = 1
        generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)
    else:  
        generator = Floodlight(surface, spectrum, math.cos(args.zenithmax), math.cos(args.zenithmin))

    tray.Add(GenerateBundles,"muons", Generator=generator,
             RunNumber=args.runid, NEvents=args.nevents,
             GCDFile=args.gcdfile,
             FromTime=dataclasses.I3Time(55380),
             ToTime=dataclasses.I3Time(55380))

    tray.AddModule('I3PropagatorModule', 'propagator',
                    PropagatorServices=make_propagators(args.rngseed),
                    RandomService=randomService, RNGStateName="RNGState")
    
    if args.useCRModel:
        tray.AddModule('I3MuonGun::WeightCalculatorModule', 'MuonWeight', Model=model, Generator=generator) 


else:
    ############nugen
    tray.AddModule("I3InfiniteSource", "source",
        prefix = args.gcdfile, 
        stream = icetray.I3Frame.DAQ
    ) 

    time = dataclasses.I3Time()
    time.set_mod_julian_time(int(args.mjd), 0, 0)

    tray.AddModule("I3MCEventHeaderGenerator", "time-gen",
        Year=time.utc_year,
        DAQTime=time.utc_daq_time,
        RunNumber=args.runid,
        IncrementEventID=True,
    )
     
    tray.AddModule(DAQCounter, "counter3",
        nevents=int(args.nevents),
    )

    try:
        tray.AddSegment(GenerateNeutrinos, 'generator',
            RandomService=randomService,
            NumEvents=args.nevents,
            SimMode=args.simmode,
            VTXGenMode=args.vtxgenmode,
            InjectionMode=args.injectionmode,
            CylinderParams=args.cylinderparams,
            AutoExtendMuonVolume=args.autoextendmuonvolume,
            Flavor=args.nuflavor,
            GammaIndex=args.gamma,
            FromEnergy = args.fromenergy,
            ToEnergy = args.toenergy,
            ZenithRange=[args.zenithmin,args.zenithmax],
            AzimuthRange=[args.azimuthmin,args.azimuthmax],
            CrossSections=args.crosssections,
            CrossSectionsPath=args.crosssectionspath,
            ParamsMap=args.paramsmap
        )
    except TypeError:
        from math import log10
        elogmin = log10(args.fromenergy)
        elogmax = log10(args.toenergy)
        
        tray.AddSegment(GenerateNeutrinos, 'generator',
            RandomService=randomService,
            NumEvents=args.nevents,
            SimMode=args.simmode,
            InjectionMode=args.injectionmode,
            CylinderParams=args.cylinderparams,
            AutoExtendMuonVolume=args.autoextendmuonvolume,
            Flavor=args.nuflavor,
            GammaIndex=args.gamma,
            EnergyLogRange=[elogmin,elogmax],
            ZenithRange=[args.zenithmin,args.zenithmax],
            AzimuthRange=[args.azimuthmin,args.azimuthmax],
            CrossSections=args.crosssections,
            CrossSectionsPath=args.crosssectionspath,
            ParamsMap=args.paramsmap
        )
        

    tray.AddSegment(PropagateMuons, 'propagator',
                    RandomService=randomService,
    ) 

    tray.AddModule(BasicCounter, "count_g",
                   Streams=[icetray.I3Frame.DAQ],
                   name="Generated Events",
    )


if len(args.sensors) == 0:
    tray.AddSegment(CLSimTraySegment,"photons",
                    GCDFile = args.gcdfile,
                    HybridMode = args.hybrid,
                    IgnoreMuons = args.hybrid,
                    gpu              = -1,
                    DisableTilt      = args.hybrid,
            )
else:
    sensors = [getattr(I3OMGeo.OMType, _) for _ in args.sensors]
    tray.AddSegment(MakePhotonsMultiSensor,
                    GCDFile=args.gcd1size,
                    Sensors=sensors,
                    UseGPUs=args.gpu,
                    EfficiencyScale=max(args.effs))
    # make pDOM hits at maximum scale as a sanity check. if the QE envelope
    # is wrong, this will raise an error
    tray.AddSegment(MakePEFromPhotons,
                    GCDFile=args.gcdfile,
                    Sensors=sensors,
                    EfficiencyScales=args.effs)
    tray.Add("Delete", Keys=['I3MCPESeriesMap'])


from icecube.gen2_sim.segments.DetectorSim import DetectorSim

inputmaps=[]
for i, Sensor in enumerate(sensors):
    SensorName = ModuleToString(Sensor)
    mapname = "I3MCPESeriesMap"+SensorName
    inputmaps.append(mapname)

tray.AddSegment(DetectorSim,
    gcdfile=args.gcdfile,
    InputMCPESeriesMaps=inputmaps,
    Sensors=sensors,
    RunID=args.runid,
    KeepMCHits=args.keepmchits,
    FilterMode=args.filtermode,
    AreaScale=args.areascale,
    VuvuzelaStartTime=args.vuvuzelastarttime,
    VuvuzelaEndTime=args.vuvuzelaendtime)

tray.Add(Calibration,"Calibration")

tray.Add(BaseProc,"BaseProc", GCDFile=args.gcdfile)#,    
#          trigger_split_config = {"IC86_SMT8": {"trigger_ids":[999], "domset_id":999},
#          "GEN2":{"trigger_ids":[999, 1003,2018], "domset_id":1000},
#          })

tray.Add(BaseReco,"BaseReco")

tray.Add("Delete", Keys=[
    'CalibratedWaveforms',
    'InIceRawData',
    'CleanInIceRawData',
    'OfflineCleanInIceRawData',
    ])

tray.AddModule("I3Writer", "writer",filename=args.outputfile)

tray.AddModule("TrashCan","trashcan")
# Execute the Tray
tray.Execute()
tray.Finish()

