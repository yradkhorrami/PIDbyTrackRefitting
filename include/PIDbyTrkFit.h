#ifndef MyPIDbyTrkFit_h
#define MyPIDbyTrkFit_h 1

#include <marlin/Processor.h>
#include <marlin/Global.h>
#include <marlinutil/HelixClass.h>
#include <marlinutil/GeometryUtil.h>
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/LCIterator.h"
#include "UTIL/Operators.h"
#include <UTIL/BitField64.h>
#include "UTIL/LCTrackerConf.h"
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>

#include "EVENT/LCStrVec.h"
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "IMPL/LCCollectionVec.h"
#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <math.h>
//#include "CalibrationHelper.h"

#include <set>
#include <vector>

class TFile;
class TH1F;
class TTree;

using namespace lcio ;
using namespace marlin ;

