// =====================================================================================
//
//       Filename:  RootTauTauEvent.h
//
//    Description:  
//
//        Version:  1.0
//        Created:  27.10.2015 16:52:08
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Ivan B. Nikolaev (ekherit), I.B.Nikolaev@inp.nsk.su
//   Organization:  Budker Institute of Nuclear Physics
//
// =====================================================================================

#pragma once

#include <vector>

#include "RootEvent/RootTuple.h"
#include "RootEvent/RootTrack.h"

//#include "RootEvent/RootPid.h"
#include "RootEvent/RootMcTruth.h"

#include "CLHEP/Vector/LorentzVector.h"

// =====================================================================================
//        Class:  RootEvent
//  Description:  Main event information supposed to used for selection
// =====================================================================================
class RootEvent : public RootTuple
{
	public:
		virtual ~RootEvent(void);
		NTuple::Item<long>    run; //run number
		NTuple::Item<long>    event; //event number 
		NTuple::Item<long>    time; //time of the event
		NTuple::Item<long>    ngood_charged_track;     //number of good charged tracks in event
		NTuple::Item<long>    ngood_neutral_track;     //number of good neutral tracks in event
		NTuple::Item<long>    ntrack;  //number of tracks (must be 2)
		RootTracks T; //track information (momentum, vertex, muon depth...)
		//RootPid Pid; //particle id for track
		RootMcTruth McTruth;
		//NTuple::Item<double>  ptsum;
		//NTuple::Item<double>  ptem;
		//NTuple::Item<double>  acop;
		//NTuple::Item<double>  acol;
		//NTuple::Item<double>  M2;

		//void init_tuple(Algorithm * algo, const char * dir, const char * title)
		//{
		//  RootTuple::init_tuple(algo,dir,title);
		//}
		virtual void bind_tuple(void)
		{
			tuple->addItem ("run", run);
			tuple->addItem ("event", event);
			tuple->addItem ("time", time);
			tuple->addItem ("channel", channel);
			tuple->addItem ("ntrack", ntrack, 0,3);
			tuple->addItem ("nctrack", nctrack, 0,10);
			tuple->addItem ("nntrack", nntrack, 0,10);
			T.add_to_tuple(tuple,ntrack); 
			McTruth.add_to_tuple(tuple,ntrack);
		};
		virtual void init(void) {};
		virtual void fill(int i,  EvtRecTrack * track);
};
