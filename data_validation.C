/*
 
  Usage:
  $ root -b -q data_validation.C(+)

  To run in compiled mode, use the "+".
 
 */

#include "rootstart.h"

#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TList.h"

#include <iostream>
#include <fstream>
#include <vector>

#define N_CHANNELS 38
#define PRESCALE 33
#define TPC_RMAX 17.44 // TPC radius (cold)
//#define TPC_FULL_MASS 43.9 //kg; used in depletion factor plots
//#define TPC_CORE_MASS 4.1  //kg; used in depletion factor plots
#define FULL_VOL_TDRIFT_MIN 40    //[us]
#define FULL_VOL_TDRIFT_MAX 310   //[us]
#define CORE_TDRIFT_MIN 135       //[us]
#define CORE_TDRIFT_MAX 235       //[us]
#define CORE_R_MAX 10             //[cm]
#define CORE2_TDRIFT_MIN 40
#define CORE2_TDRIFT_MAX 334.6
#define CORE2_R_MAX 14.14         // = sqrt(200)
#define VDRIFT 0.93               //[mm/us]
#define LAR_DENSITY 1.4           //[g/cm^3]

#define FV_TDRIFT_MIN 40
#define FV_TDRIFT_MAX 334.6
#define S1_MIN 60.
#define S1_MAX 460.

// Open output file
TString outfile_name = "NR_for_ML.root";
TFile* outfile = new TFile(outfile_name, "RECREATE");

Bool_t load_allpulses = true;
Bool_t load_masas_xy = true;
Bool_t load_xylocator_xy = true;
Bool_t load_aww_xy = true;
Bool_t load_veto = true;

// Enum to help keep track of cuts
const Int_t nCuts = 24;
enum EventCut_t {
  NOCUTS,        NCHANNEL,     BASELINE,      LIVEINHIBIT,    LONGWAIT,                      //BASIC
  VETO_PRESENT,  VETO_PROMPT,  VETO_DELAYED,  VETO_PREPROMPT, VETO_MUON,  VETO_COSMOGENIC,   //VETO PHYSICS
  NPULSES,       TRIGGERTIME,  S1SATURATION,  S1MAXFRAC,      S2VALID,    S2MINIMUM,         //TPC PHYSICS
  S2MAXIMUM,     Z_FIDUCIAL,   S1RANGE,                                                      //TPC PHYSICS
  RADIUS16,      S2OVERS1,     RADIUS10,      FASTf90                                        //TEST CUTS
}; // enum starts from 0. ie, NOCUTS=0, NCHANNEL=1, ...

const EColor colors[26] = {kBlack, kGray, kYellow, kRed, kMagenta, kBlue, kCyan, kGreen,
                           EColor(kOrange-5), EColor(kPink-5), EColor(kViolet-5), EColor(kAzure-5), EColor(kTeal-5), EColor(kSpring-5), 
                           EColor(kYellow+2), EColor(kRed+2), EColor(kMagenta+2), EColor(kBlue+2), EColor(kCyan+2), EColor(kGreen+2),
                           EColor(kYellow-9), EColor(kRed-9), EColor(kMagenta-9), EColor(kBlue-9), EColor(kCyan-9), EColor(kGreen-9)};

struct TPCEvent {
  // variables from main file
  // SLAD RunXXXXXX.root
  //events
  Int_t tpc_run_id;
  Int_t tpc_subrun_id;
  Int_t tpc_event_id;
  UInt_t tpc_dig_sum;
  //gps
  Int_t tpc_gps_coarse;
  Int_t tpc_gps_fine;
  //long_lifetime
  Double_t live_time;
  Double_t inhibit_time;
  //logbook
  //trigger
  UInt_t trigger_type;
  
  //nchannels
  Int_t nchannels;
  //baseline
  Short_t baseline_not_found;
  //npulses
  Int_t npulses;
  Int_t has_s3;
  Int_t has_s1echo;
  //tdrift
  Double_t tdrift;
  //s1
  Float_t total_s1;
  Float_t total_s1_corr;
  Float_t total_s1_top;
  Float_t total_s1_bottom;
  Float_t total_s1_long_us_integrals[170];
  Float_t total_s2_long_us_integrals[170]; 
  //s1_saturation
  Int_t s1_saturated;
  //s1_time
  Float_t s1_start_time;
  //s1_f90
  Float_t total_f90;
  //s1_fraction
  Int_t s1_max_chan;
  Float_t s1_max_frac;
  //max_s1_frac_cut
  Float_t max_s1_frac_cut_threshold99;
  Int_t max_s1_frac_cut_exceeds99;
  //s2
  Float_t total_s2;
  //s2_saturation
  Int_t s2_saturated;
  //s2_f90
  Float_t total_s2_f90_fixed;
  //bary_s2
  
  // SLAD RunXXXXXX_s2.root
  //s2_fraction
  //Int_t s2_max_chan;
  //Float_t s2_max_frac;
  //Float_t s2_ch_frac[N_CHANNELS];
  
  // SLAD RunXXXXXX_masas_xy.root
  // variables from xy file
  Float_t x_masas;
  Float_t y_masas;
  Float_t xycorr_factor_masas;
  Float_t r_masas;
  Float_t theta_masas;

  Float_t pulse_x_masas[100];
  Float_t pulse_y_masas[100];
  Float_t pulse_xycorr_factor_masas[100];
  Float_t pulse_r_masas[100];
  Float_t pulse_theta_masas[100];
  
  // SLAD RunXXXXXX_xylocator_xy.root
  // variables from xy file
  Float_t x_jasons;
  Float_t y_jasons;
  Float_t xycorr_factor_jasons;
  Float_t r_jasons;
  Float_t theta_jasons;

  Float_t pulse_x_jasons[100];
  Float_t pulse_y_jasons[100];
  Float_t pulse_xycorr_factor_jasons[100];
  Float_t pulse_r_jasons[100];
  Float_t pulse_theta_jasons[100];
 
  // SLAD RunXXXXXX_xylocator_xy.root
  // variables from xy file
  Float_t x_aww;
  Float_t y_aww;
  Float_t xycorr_factor_aww;
  Float_t r_aww;
  Float_t theta_aww;
  
  Float_t pulse_x_aww[100];
  Float_t pulse_y_aww[100];
  Float_t pulse_xycorr_factor_aww[100];
  Float_t pulse_r_aww[100];
  Float_t pulse_theta_aww[100];
  
  // SLAD RunXXXXXX_allpulses.root
  //pulse_info
  Float_t pulse_start_time[100];
  Float_t pulse_max_ch[100];
  Float_t pulse_us_integrals[7*100];
  Float_t pulse_peak_time[100];
  Float_t pulse_peak_amp[100];
  
  // SLAD RunXXXXXX_veto.root
  Int_t veto_run_id;
  Int_t veto_event_id;
  Int_t veto_present;
  Float_t veto_lsv_total_charge;
  Float_t veto_wt_total_charge;
  Int_t	veto_nclusters; 
  std::vector<float>* veto_roi_lsv_charge_vec;
  std::vector<float>* veto_slider_lsv_charge_vec;
  std::vector<float>* veto_cluster_dtprompt_vec;
  std::vector<float>* veto_cluster_charge_vec;

  // Variables to be calculated on the fly
  Float_t s1_prompt;
  Float_t s1_late;
  Float_t s1_tbacorr;
  Float_t s1_prompt_tbacorr;
  Float_t s1_late_tbacorr;
  Float_t x, y, z, r, xycorr_factor;
  Float_t total_s2_xycorr;
  Float_t s1mf_threshold;

  Float_t nr_median_s2;
  Float_t nr_rms_s2;
  Float_t nr_sigma_s2;

  Float_t nr_s2overs1;
  Float_t s2overs1;

  Float_t kr85mpeaktime;
  Float_t p0f5000;

  Float_t muon_dt;   // Time since previous muon
};

//------------------------------------------------------------------------------
// Set the branch address for the variables in TPCEvent.
// To increase speed, we disable all branches by default, and turn on only the
// ones we want to use.

void load_tpctree(TChain* tpc_events, TPCEvent& e) {
  tpc_events->SetBranchStatus("*",0); //disable all
  
  // SLAD RunXXXXXX.root
  // events
  tpc_events->SetBranchStatus("events.run_id", 1);
  tpc_events->SetBranchAddress("run_id", &e.tpc_run_id);

  tpc_events->SetBranchStatus("events.subrun_id", 1);
  tpc_events->SetBranchAddress("subrun_id", &e.tpc_subrun_id);

  tpc_events->SetBranchStatus("events.event_id", 1);
  tpc_events->SetBranchAddress("event_id", &e.tpc_event_id);

  tpc_events->SetBranchStatus("events.tpc_digital_sum", 1);
  tpc_events->SetBranchAddress("tpc_digital_sum", &e.tpc_dig_sum);
  
  // gps
  tpc_events->SetBranchStatus("gps.gps_coarse", 1);
  tpc_events->SetBranchAddress("gps.gps_coarse", &e.tpc_gps_coarse);
  
  tpc_events->SetBranchStatus("gps.gps_fine", 1);
  tpc_events->SetBranchAddress("gps.gps_fine", &e.tpc_gps_fine);
  
  // long_lifetime
  tpc_events->SetBranchStatus("long_lifetime.lifetime",1);
  tpc_events->SetBranchAddress("long_lifetime.lifetime", &e.live_time);

  tpc_events->SetBranchStatus("long_lifetime.inhibittime",1);
  tpc_events->SetBranchAddress("long_lifetime.inhibittime", &e.inhibit_time);
  
  //logbook
  //trigger_type
  tpc_events->SetBranchStatus("trigger_type", 1);
  tpc_events->SetBranchAddress("trigger_type", &e.trigger_type);
  
  //nchannels
  tpc_events->SetBranchStatus("nchannel.nchannel", 1);
  tpc_events->SetBranchAddress("nchannel.nchannel", &e.nchannels);

  // baseline
  tpc_events->SetBranchStatus("baseline.SumChannelHasNoBaseline",1);
  tpc_events->SetBranchAddress("baseline.SumChannelHasNoBaseline", &e.baseline_not_found);

  // npulses
  tpc_events->SetBranchStatus("npulses.n_phys_pulses",1);
  tpc_events->SetBranchAddress("npulses.n_phys_pulses", &e.npulses);

  tpc_events->SetBranchStatus("npulses.has_s3",1);
  tpc_events->SetBranchAddress("npulses.has_s3", &e.has_s3);

  tpc_events->SetBranchStatus("npulses.has_s1echo",1);
  tpc_events->SetBranchAddress("npulses.has_s1echo", &e.has_s1echo);
  
  // tdrift
  tpc_events->SetBranchStatus("tdrift.tdrift", 1);
  tpc_events->SetBranchAddress("tdrift.tdrift", &e.tdrift);

  // s1
  tpc_events->SetBranchStatus("s1.total_s1", 1);
  tpc_events->SetBranchAddress("s1.total_s1", &e.total_s1);
  
  tpc_events->SetBranchStatus("s1.total_s1_corr", 1);
  tpc_events->SetBranchAddress("s1.total_s1_corr", &e.total_s1_corr);
  
  tpc_events->SetBranchStatus("s1.total_s1_top", 1);
  tpc_events->SetBranchAddress("s1.total_s1_top", &e.total_s1_top);
  
  tpc_events->SetBranchStatus("s1.total_s1_bottom", 1);
  tpc_events->SetBranchAddress("s1.total_s1_bottom", &e.total_s1_bottom);

  tpc_events->SetBranchStatus("s1.total_s1_long_us_integrals", 1);
  tpc_events->SetBranchAddress("s1.total_s1_long_us_integrals", &e.total_s1_long_us_integrals);
  
  // s1_saturation
  tpc_events->SetBranchStatus("s1_saturation.is_saturated_pulse0", 1);
  tpc_events->SetBranchAddress("s1_saturation.is_saturated_pulse0", &e.s1_saturated);

  // s1_time
  tpc_events->SetBranchStatus("s1_time.s1_start_time", 1);
  tpc_events->SetBranchAddress("s1_time.s1_start_time", &e.s1_start_time);

  // s1_f90
  tpc_events->SetBranchStatus("s1_f90.total_f90", 1);
  tpc_events->SetBranchAddress("s1_f90.total_f90", &e.total_f90);
  
  // s1_fraction
  tpc_events->SetBranchStatus("s1_fraction.s1_max_chan", 1);
  tpc_events->SetBranchAddress("s1_fraction.s1_max_chan", &e.s1_max_chan);
  
//  tpc_events->SetBranchStatus("s1_fraction.s1_max_frac", 1);
//  tpc_events->SetBranchAddress("s1_fraction.s1_max_frac", &e.s1_max_frac);

  tpc_events->SetBranchStatus("s1_fraction.s1_prompt_max_frac", 1);
  tpc_events->SetBranchAddress("s1_fraction.s1_prompt_max_frac", &e.s1_max_frac);
  
  // max_s1_frac_cut
  tpc_events->SetBranchStatus("max_s1_frac_cut.max_s1_frac_cut_threshold99", 1);
  tpc_events->SetBranchAddress("max_s1_frac_cut.max_s1_frac_cut_threshold99", &e.max_s1_frac_cut_threshold99);
  
  tpc_events->SetBranchStatus("max_s1_frac_cut.max_s1_frac_cut_exceeds99", 1);
  tpc_events->SetBranchAddress("max_s1_frac_cut.max_s1_frac_cut_exceeds99", &e.max_s1_frac_cut_exceeds99);
  
  // s2
  tpc_events->SetBranchStatus("s2.total_s2", 1);
  tpc_events->SetBranchAddress("s2.total_s2", &e.total_s2);
  
  tpc_events->SetBranchStatus("s2.total_s2_long_us_integrals", 1);
  tpc_events->SetBranchAddress("s2.total_s2_long_us_integrals", &e.total_s2_long_us_integrals);


  // s2_saturation
  tpc_events->SetBranchStatus("s2_saturation.is_saturated_pulse1", 1);
  tpc_events->SetBranchAddress("s2_saturation.is_saturated_pulse1", &e.s2_saturated);


  // s2_f90
  tpc_events->SetBranchStatus("s2_f90.total_s2_f90", 1);
  tpc_events->SetBranchAddress("s2_f90.total_s2_f90", &e.total_s2_f90_fixed);

  if (load_masas_xy) {
    // SLAD RunXXXXXX_masas_xy.root
    // masa_x
    tpc_events->SetBranchStatus("masas_x", 1);
    tpc_events->SetBranchAddress("masas_x", &e.x_masas);
    
    // masa_y
    tpc_events->SetBranchStatus("masas_y", 1);
    tpc_events->SetBranchAddress("masas_y", &e.y_masas);
    
    // masas_chi2
    
    // masas_xycorr_factor
    tpc_events->SetBranchStatus("masas_xycorr_factor", 1);
    tpc_events->SetBranchAddress("masas_xycorr_factor", &e.xycorr_factor_masas);
    
    // masa_r
    tpc_events->SetBranchStatus("masas_r", 1);
    tpc_events->SetBranchAddress("masas_r", &e.r_masas);
    
    // masa_theta
    tpc_events->SetBranchStatus("masas_theta", 1);
    tpc_events->SetBranchAddress("masas_theta", &e.theta_masas);
    
    // masa_x
    tpc_events->SetBranchStatus("allpulses_x", 1);
    tpc_events->SetBranchAddress("allpulses_x", &e.pulse_x_masas);
    
    // masa_y
    tpc_events->SetBranchStatus("allpulses_y", 1);
    tpc_events->SetBranchAddress("allpulses_y", &e.pulse_y_masas);
    
    // masas_chi2
    
    // masas_xycorr_factor
    tpc_events->SetBranchStatus("allpulses_xycorr_factor", 1);
    tpc_events->SetBranchAddress("allpulses_xycorr_factor", &e.pulse_xycorr_factor_masas);
    
    // masa_r
    tpc_events->SetBranchStatus("allpulses_r", 1);
    tpc_events->SetBranchAddress("allpulses_r", &e.pulse_r_masas);
    
    // masa_theta
    tpc_events->SetBranchStatus("allpulses_theta", 1);
    tpc_events->SetBranchAddress("allpulses_theta", &e.pulse_theta_masas);
  }

  if (load_xylocator_xy) {
    // SLAD RunXXXXXX_xylocator_xy.root
    // xyl_SCM
    
    // xyl_best_x
    tpc_events->SetBranchStatus("xyl_best_x", 1);
    tpc_events->SetBranchAddress("xyl_best_x", &e.x_jasons);
    
    // xyl_best_y
    tpc_events->SetBranchStatus("xyl_best_y", 1);
    tpc_events->SetBranchAddress("xyl_best_y", &e.y_jasons);
    
    // xyl_best_chi2
    
    // xylocator xy correction factor
    tpc_events->SetBranchStatus("xyl_best_xycorr_factor", 1);
    tpc_events->SetBranchAddress("xyl_best_xycorr_factor", &e.xycorr_factor_jasons);
    
    // xyl_best_r
    tpc_events->SetBranchStatus("xyl_best_r", 1);
    tpc_events->SetBranchAddress("xyl_best_r", &e.r_jasons);
    
    // xyl_best_theta
    tpc_events->SetBranchStatus("xyl_best_theta", 1);
    tpc_events->SetBranchAddress("xyl_best_theta", &e.theta_jasons);
    
    
    // xyl_best_x
    tpc_events->SetBranchStatus("allpulses_xyl_x", 1);
    tpc_events->SetBranchAddress("allpulses_xyl_x", &e.pulse_x_jasons);
    
    // xyl_best_y
    tpc_events->SetBranchStatus("allpulses_xyl_y", 1);
    tpc_events->SetBranchAddress("allpulses_xyl_y", &e.pulse_y_jasons);
    
    // xyl_best_chi2
    
    // xylocator xy correction factor
    tpc_events->SetBranchStatus("allpulses_xyl_xycorr_factor", 1);
    tpc_events->SetBranchAddress("allpulses_xyl_xycorr_factor", &e.pulse_xycorr_factor_jasons);
    
    // xyl_best_r
    tpc_events->SetBranchStatus("allpulses_xyl_r", 1);
    tpc_events->SetBranchAddress("allpulses_xyl_r", &e.pulse_r_jasons);
    
    // xyl_best_theta
    tpc_events->SetBranchStatus("allpulses_xyl_theta", 1);
    tpc_events->SetBranchAddress("allpulses_xyl_theta", &e.pulse_theta_jasons);
  }

  if (load_aww_xy) {
    // SLAD RunXXXXXX_aww_xy.root
    // aww_x
    tpc_events->SetBranchStatus("aww_x", 1);
    tpc_events->SetBranchAddress("aww_x", &e.x_aww);
    
    // aww_y
    tpc_events->SetBranchStatus("aww_y", 1);
    tpc_events->SetBranchAddress("aww_y", &e.y_aww);
    
    // aww_chi2
    
    // aww_xycorr_factor
    tpc_events->SetBranchStatus("aww_xycorr_factor", 1);
    tpc_events->SetBranchAddress("aww_xycorr_factor", &e.xycorr_factor_aww);
    
    // aww_r
    tpc_events->SetBranchStatus("aww_r", 1);
    tpc_events->SetBranchAddress("aww_r", &e.r_aww);
    
    // aww_theta
    tpc_events->SetBranchStatus("aww_theta", 1);
    tpc_events->SetBranchAddress("aww_theta", &e.theta_aww);
    
    // aww_x
    tpc_events->SetBranchStatus("allpulses_x", 1);
    tpc_events->SetBranchAddress("allpulses_x", &e.pulse_x_aww);
    
    // aww_y
    tpc_events->SetBranchStatus("allpulses_y", 1);
    tpc_events->SetBranchAddress("allpulses_y", &e.pulse_y_aww);
    
    // aww_chi2
    
    // aww_xycorr_factor
    tpc_events->SetBranchStatus("allpulses_xycorr_factor", 1);
    tpc_events->SetBranchAddress("allpulses_xycorr_factor", &e.pulse_xycorr_factor_aww);
    
    // aww_r
    tpc_events->SetBranchStatus("allpulses_r", 1);
    tpc_events->SetBranchAddress("allpulses_r", &e.pulse_r_aww);
    
    // aww_theta
    tpc_events->SetBranchStatus("allpulses_theta", 1);
    tpc_events->SetBranchAddress("allpulses_theta", &e.pulse_theta_aww);
  }

  if (load_allpulses) {
    // SLAD RunXXXXXX_allpulses.root
    // pulse_info
    tpc_events->SetBranchStatus("pulse_info.pulse_info_start_time", 1);
    tpc_events->SetBranchAddress("pulse_info.pulse_info_start_time", e.pulse_start_time);
    
    tpc_events->SetBranchStatus("pulse_info.pulse_info_max_chan", 1);
    tpc_events->SetBranchAddress("pulse_info.pulse_info_max_chan", e.pulse_max_ch);

    tpc_events->SetBranchStatus("pulse_info.pulse_info_us_integrals", 1);
    tpc_events->SetBranchAddress("pulse_info.pulse_info_us_integrals", e.pulse_us_integrals);

    tpc_events->SetBranchStatus("pulse_info.pulse_info_peak_time", 1);
    tpc_events->SetBranchAddress("pulse_info.pulse_info_peak_time", e.pulse_peak_time);
    
    tpc_events->SetBranchStatus("pulse_info.pulse_info_peak_amp", 1);
    tpc_events->SetBranchAddress("pulse_info.pulse_info_peak_amp", e.pulse_peak_amp);
  }

  if (load_veto) {
    // SLAD RunXXXXXX_veto.root
    tpc_events->SetBranchStatus("veto_run_id", 1);
    tpc_events->SetBranchAddress("veto_run_id", &e.veto_run_id);

    tpc_events->SetBranchStatus("veto_nclusters", 1);
    tpc_events->SetBranchAddress("veto_nclusters", &e.veto_nclusters);    

    tpc_events->SetBranchStatus("veto_event_id", 1);
    tpc_events->SetBranchAddress("veto_event_id", &e.veto_event_id);
    
    tpc_events->SetBranchStatus("veto_present", 1);
    tpc_events->SetBranchAddress("veto_present", &e.veto_present);
    
    tpc_events->SetBranchStatus("veto_lsv_total_charge", 1);
    tpc_events->SetBranchAddress("veto_lsv_total_charge", &e.veto_lsv_total_charge);
    
    tpc_events->SetBranchStatus("veto_wt_total_charge", 1);
    tpc_events->SetBranchAddress("veto_wt_total_charge", &e.veto_wt_total_charge);
  
    /*
      ROIs
      new
      0  | [-50  ,250]   ns prompt
      1  | [-2000,-1700] ns pre-prompt bg
      2  | [40000,40300] ns post-prompt bg
      old
      3  | [-10  ,200]   ns prompt
      4  | [-2000,-1790] ns pre-prompt bg
      5  | [40000,40210] ns post-prompt bg
   
      SLIDERs
      new
      0  | [0 ns ,end]  500 ns post-prompt
      1  | [start,0 ns] 500 ns pre-prompt
      old
      2  | [0      ,8800 ns] 300 ns early post-prompt
      3  | [8800 ns,end]     300 ns early pre-prompt
    */
  
    e.veto_roi_lsv_charge_vec = 0;
   e.veto_slider_lsv_charge_vec = 0;
    e.veto_cluster_dtprompt_vec = 0;
    e.veto_cluster_charge_vec = 0;
  
    tpc_events->SetBranchStatus("veto_roi_lsv_charge_vec", 1);
    tpc_events->SetBranchAddress("veto_roi_lsv_charge_vec", &e.veto_roi_lsv_charge_vec);
  
    tpc_events->SetBranchStatus("veto_slider_lsv_charge_vec", 1);
    tpc_events->SetBranchAddress("veto_slider_lsv_charge_vec", &e.veto_slider_lsv_charge_vec);

	
    tpc_events->SetBranchStatus("veto_cluster_dtprompt_vec", 1);
    tpc_events->SetBranchAddress("veto_cluster_dtprompt_vec", &e.veto_cluster_dtprompt_vec);
  
    tpc_events->SetBranchStatus("veto_cluster_charge_vec", 1);
    tpc_events->SetBranchAddress("veto_cluster_charge_vec", &e.veto_cluster_charge_vec);
  }
} //load_tpctree

Float_t s1_tbacorr(Float_t s1, Float_t tbasym){
  
  Float_t s1_par0 = 0.948485;
  Float_t s1_par1 = -0.793215;
  Float_t s1_par2 = -2.26064;
  Float_t s1_par3 = 6.84532;
  Float_t s1_par4 = 73.8491;
  Float_t s1_par5 = 127.36;
  Float_t s1_corr = 0.;

  if (tbasym <= -0.35)   s1_corr = s1/(s1_par0 + s1_par1*-0.35 + s1_par2*pow(-0.35,2) + s1_par3*pow(-0.35,3) + s1_par4*pow(-0.35,4) + s1_par5*pow(-0.35,5));
  else if(tbasym < 0.15) s1_corr = s1/(s1_par0 + s1_par1*tbasym + s1_par2*pow(tbasym,2) + s1_par3*pow(tbasym,3) + s1_par4*pow(tbasym,4) + s1_par5*pow(tbasym,5));
  else                   s1_corr = s1/(s1_par0 + s1_par1*0.15 + s1_par2*pow(0.15,2) + s1_par3*pow(0.15,3) + s1_par4*pow(0.15,4) + s1_par5*pow(0.15,5));

  return s1_corr;
}

//---------------------------------------
// AAr G2 trigger has pre-scale thresholds, and the thresholds are run dependent.
//
// The function below check to see if an events is pre-scaled
// The pre-scale factor is 33, meaning only 1/33 of AAr events is saved
// If this event is pre-scaled, you can fill the same event 33 times
//
// Note, below check is only valid for 50d AAr (calibration run are more complicated)
//----------------------------------------


Bool_t isPrescaled(TPCEvent const& e) {
  
  UInt_t thres[2] = {999999, 999999};
  if (e.tpc_run_id >= 5372 && e.tpc_run_id <= 8433) {
    if (e.tpc_run_id < 7854) {
      thres[0] = 999999;
      thres[1] = 999999;
    }
    else if (e.tpc_run_id >= 7854 && e.tpc_run_id < 8276){
      thres[0] = 360;
      thres[1] = 999999;
    }
    else{
      thres[0] = 360;
      thres[1] = 1500;
    }
  }
  else return false;
  
  if (e.tpc_dig_sum >= thres[0] && e.tpc_dig_sum < thres[1]) return true;
  else return false;
}



//------------------------------------------------------------------------------
// Define cuts.

//--- Quality cuts
Bool_t cx_nchan(TPCEvent const& e)        { return e.nchannels == N_CHANNELS; }
Bool_t cx_baseline(TPCEvent const& e)     { return e.baseline_not_found == false; }
Bool_t cx_old_event_dt(TPCEvent const& e) { return (e.live_time + e.inhibit_time) > 1.35E-3; } // used in analyses up to 70d.
Bool_t cx_event_dt(TPCEvent const& e)     { return e.live_time > 400.E-6; }
Bool_t cx_file_io(TPCEvent const& e)      { return e.live_time < 1.; } // only for AAr
Bool_t cx_veto_present(TPCEvent const& e) { return e.veto_present > 0; }

//--- Old veto cuts (for AAr)
Bool_t cx_old_veto_prompt(TPCEvent const& e)  { return e.veto_present && e.veto_roi_lsv_charge_vec->at(3) < 10.; }
Bool_t cx_old_veto_delayed(TPCEvent const& e) { return
    e.veto_present && e.veto_slider_lsv_charge_vec->at(2) < 80. && e.veto_slider_lsv_charge_vec->at(3) < 110. && e.veto_wt_total_charge < 200.; }
    
//--- Veto cuts for UAr
Bool_t cx_veto_prompt(TPCEvent const& e)     { return e.veto_present && !(e.veto_roi_lsv_charge_vec->at(0) > 1.); }
Bool_t cx_veto_delayed(TPCEvent const& e)    { return e.veto_present &&
    !((e.veto_run_id < 12638 && e.veto_slider_lsv_charge_vec->at(0) > 3.) || (e.veto_run_id >= 12638 && e.veto_slider_lsv_charge_vec->at(0) > 6.)); }
Bool_t cx_veto_preprompt(TPCEvent const& e)  { return e.veto_present && !(e.veto_slider_lsv_charge_vec->at(1) > 3.); }
Bool_t cx_veto_muon(TPCEvent const& e)       { return e.veto_present && !(e.veto_lsv_total_charge > 2000. || e.veto_wt_total_charge > 400.); }
Bool_t cx_veto_cosmogenic(TPCEvent const& e) { return e.veto_present && e.muon_dt > 2.; }

//--- For estimating Veto prompt cut acceptance
Bool_t cx_veto_prompt_acc(TPCEvent const& e)     { return e.veto_present && e.veto_roi_lsv_charge_vec->at(1) > 1.; }
Bool_t cx_veto_prompt_acc_sys(TPCEvent const& e) { return e.veto_present && e.veto_roi_lsv_charge_vec->at(2) > 1.; }

//--- Physics cuts
Bool_t cx_single_scatter(TPCEvent const& e)  { return e.npulses == 2 || (e.npulses == 3 && e.has_s3); }
Bool_t cx_trg_time(TPCEvent const& e) { return
  ((e.tpc_run_id >= -999 && e.tpc_run_id < 7344 && e.s1_start_time >= -0.25 && e.s1_start_time <= -0.15) ||
   (e.tpc_run_id >= 7344 && e.tpc_run_id < 7641 && e.s1_start_time >= -4.10 && e.s1_start_time <= -4.00) ||
   (e.tpc_run_id >= 7641 && e.tpc_run_id < 999999 && e.s1_start_time >= -6.10 && e.s1_start_time <= -6.00)); }
Bool_t cx_s1_sat(TPCEvent const& e)   { return e.s1_saturated == 0; }
Bool_t cx_s1_mf(TPCEvent const& e)    { return e.s1_max_frac < e.s1mf_threshold; }
Bool_t cx_s2_f90(TPCEvent const& e)   { return e.total_s2_f90_fixed < 0.2; }
Bool_t cx_s2_size(TPCEvent const& e)  { return e.r > 0. && e.r < 20. && e.total_s2_xycorr > 100.; }
Bool_t cx_max_s2(TPCEvent const& e)   { return e.r > 0. && e.r < 20. && e.total_s2_xycorr < 8000.; } // NEW
Bool_t cx_s1_range(TPCEvent const& e) { return e.total_s1_corr > S1_MIN && e.total_s1_corr < S1_MAX; }
Bool_t cx_fiducial(TPCEvent const& e) { return e.tdrift > FV_TDRIFT_MIN && e.tdrift < FV_TDRIFT_MAX; }
Bool_t cx_r(TPCEvent const& e, Float_t r_cx) { return e.r > 0 && e.r < r_cx; }
//Bool_t cx_s2overs1(TPCEvent const& e) { return e.nr_sigma_s2 < -0.1; }
Bool_t cx_s2overs1(TPCEvent const& e) { return e.s2overs1 < e.nr_s2overs1; }
Bool_t cx_f90_fast(TPCEvent const& e) { return e.total_f90 > 0.05; }  //2015-09-05 AFan: Changed from .15 to .05 for null field runs. wasn't really being used for field on analysis anyway

// Select core of TPC
Bool_t cx_core(TPCEvent const& e)     { return e.tdrift > CORE_TDRIFT_MIN && e.tdrift < CORE_TDRIFT_MAX && e.r >= 0 && e.r < CORE_R_MAX; }
Bool_t cx_core2(TPCEvent const& e)    { return e.tdrift > CORE2_TDRIFT_MIN && e.tdrift < CORE2_TDRIFT_MAX && e.r >= 0 && e.r < CORE2_R_MAX; }

// Use maximum TPC volume for good events. Outside this, weird reconstruction effects may appear.
Bool_t cx_full_volume(TPCEvent const& e) { return e.tdrift > FULL_VOL_TDRIFT_MIN && e.tdrift < FULL_VOL_TDRIFT_MAX && e.r > 0 && e.r < TPC_RMAX; }

// Cuts for Kr85m search
Bool_t cx_kr85m_1p(TPCEvent const& e)        { return (e.npulses >= 2 && e.pulse_start_time[1] - e.pulse_start_time[0] > 5.); }
Bool_t cx_kr85m_2p(TPCEvent const& e)        { return (e.npulses >= 3 && e.pulse_start_time[1] - e.pulse_start_time[0] <= 5. &&
                                                       e.pulse_start_time[2] - e.pulse_start_time[0] > 5.); }
Bool_t cx_kr85m_peaktime(TPCEvent const& e)  { return e.kr85mpeaktime > 0.05 && e.kr85mpeaktime < 5; }
Bool_t cx_kr85m_f5000(TPCEvent const& e)     { return e.p0f5000 > 0.90; }
Bool_t cx_kr85m_box(TPCEvent const& e)       { return e.npulses >= 2 && e.total_s1 > 400. && e.total_s1 < 10.E+3 && e.total_f90 < 0.2; }

Bool_t cx_f90_s2s1_s1range(TPCEvent const& e){ return e.total_s1_corr > 200 && e.total_s1_corr < 400; }

//------------------------------------------------------------------------------
// Main event loop is contained here.
void event_loop(TChain* tpc_events, TString type = "") {
  TH1::SetDefaultSumw2(kTRUE);
  
  // Set the branch addresses.
  TPCEvent e;
  load_tpctree(tpc_events, e);

  Bool_t aar  = (type == "aar");
  Bool_t uar  = (type == "uar");
  Bool_t ambe = (type == "ambe" || type == "ambebg");
	Float_t new_v;
  TTree *t2 = new TTree("t2","data from histogram");  
  t2->Branch("point",&new_v,"s1 form");

// Load max_s1_frac threshold file so can apply other versions of S1 max frac cut on the fly.
  TFile* s1mf_file = new TFile("max_frac_cut_fixed_acceptance_full_stats.root");
  if(s1mf_file->IsZombie() || !s1mf_file) std::cout << "Bad S1mf_file" << std::endl;
  TH2F* h_s1mf_thresholds = (TH2F*) s1mf_file->Get("s1pmf_c95"); // Load 95%ile cut.
  
  // Load 1-sigma contours for getting NR from AmBe data
  TFile* dmbox_file = new TFile("dmbox.root");
  if(dmbox_file->IsZombie() || !dmbox_file) std::cout << "Bad dmbox_file (getting NR 99% acc from here)" << std::endl;
  TGraph* ambe_acc99 = (TGraph*) dmbox_file->Get("NRacceptances/g_acceptance_99");

  // Define histograms
  outfile->cd();
  outfile->mkdir(type);
  outfile->cd(type);
  
  // EventCounter
  TDirectory* metadata_dir = gDirectory->mkdir("metadata");
  metadata_dir->cd();
  const TString CutLabels[nCuts] = {
    "No cuts",      "# of channels",  "Baseline",      "Live + inhibit time",  "Long wait",
    "Veto present", "LSV prompt",     "LSV delayed",   "LSV pre-prompt",       "Muon",       "Cosmogenic",
    "# of pulses",  "Trigger time",   "S1 saturation", "S1 maximum fraction",  "S2 valid",   "S2 minimum",
    "Drift time fiducialization",     "S1 range",
    "Radius 16 cm", "S2/S1",          "Radius 10 cm",  "Fast f90"
  };
  TH1I* h_EventCounter = new TH1I("h_EventCounter","Event Counter; ; Number of events", nCuts, -0.5, nCuts-0.5);
  for (Int_t ieb = 1; ieb <= nCuts; ++ieb) h_EventCounter->GetXaxis()->SetBinLabel(ieb,CutLabels[ieb-1].Data());
  
  outfile->cd(type);
  
  // Other: S2 distribution, s2_f90 distribution, tdrift, s1 max frac, etc
  TDirectory* other_dir = gDirectory->mkdir("other");
  other_dir->cd();
  TH1F* h_s1_start     = new TH1F("h_s1_start",     "; S1 start time [#mus]", 1000, -10., 10.);
  TH1F* h_tdrift       = new TH1F("h_tdrift",       "; t_{drift} [#mus]", 200, 0., 400.);
  TH1F* h_tdrift_lowS1 = new TH1F("h_tdrift_lowS1", "; t_{drift} [#mus]", 200, 0., 400.);
  TH1F* h_s1saturation = new TH1F("h_s1saturation", "; S1 [PE]", 1000, 0., 10.E+3);
  TH2F* h_xy_masas     = new TH2F("h_xy_masas",     "; x [cm]; y [cm]", 40, -20., 20., 40, -20., 20.);
  TH2F* h_xy_jasons    = new TH2F("h_xy_jasons",    "; x [cm]; y [cm]", 40, -20., 20., 40, -20., 20.);
  TH2F* h_rz_masas     = new TH2F("h_rz_masas",     "; r^{2}/r_{TPC} [cm]; -t_{drift} [-#mus]", 20, 0., 20., 100, -390., 10.);
  TH2F* h_rz_jasons    = new TH2F("h_rz_jasons",    "; r^{2}/r_{TPC} [cm]; -t_{drift} [-#mus]", 20, 0., 20., 100, -390., 10.);
  TH2F* h_rz           = new TH2F("h_rz",           "; r^{2}/r_{TPC} [cm]; -t_{drift} [-#mus]", 20, 0., 20., 100, -390., 10.);
  TH2F* h_S2xycorr     = new TH2F("h_S2xycorr",     "; S1 [PE]; S2_{xycorr} [PE]", 500, 0., 500., 1000, 0., 50.E+3);
  TH1F* h_s1_badxy     = new TH1F("h_s1_badxy",     "; S1 [PE]", 1000, 0., 10.E+3);
  TH1F* h_s2_badxy     = new TH1F("h_s2_badxy",     "; S2 [PE]", 1000, 0., 100.E+3);
  TH1F* h_tdrift_badxy = new TH1F("h_tdrift_badxy", "; t_{drift} [#mus]", 200, 0., 400.);
  TH1F* h_S2f90        = new TH1F("h_S2f90",        "; S2_f90", 1000, 0., 1.);
  TH1F* h_rdiff        = new TH1F("h_rdiff",        "; #Deltar [cm]", 400, -20., 20.);
  Double_t norm = 1;
  TH1F* h_s1_waveform  = new TH1F("h_s1_form",      ";s1 long integrals, ", 500, 0., 9000000.);
  TH1F* h_s2_waveform  = new TH1F("h_s2_form",      ";s2 long integrals, ", 500, 0., 9000000.);
  outfile->cd(type);

  // DMS
  TDirectory* dms_dir = gDirectory->mkdir("dms");
  dms_dir->cd();
  const Int_t ndms = 12;
  TH2F** h_DMS = new TH2F*[ndms];
  const Int_t nbins_s1_dms = 200, nbins_f90_dms = 100;
  const Float_t s1_min_dms = 0., s1_max_dms = 1000., f90_min_dms = 0., f90_max_dms = 1.;
  const TString h_DMSLabels[ndms] = {
    "Dark Matter Search (without veto)","Dark Matter Search",
    "Dark Matter Search (veto prompt only)", "Dark Matter Search (veto delayed only)", "Dark Matter Search (veto preprompt only)", 
    "Dark Matter Search (veto muon only)", "Dark Matter Search (veto cosmogenic only)",
    "Dark Matter Search (no veto cuts, no s1mf cut)", "Dark Matter Search (pass veto cuts, no s1mf cut)",
    "Dark Matter Search (S2/S1, deep R)", "Dark Matter Search (deep R)", "Dark Matter Search (core)"
  };
  for (Int_t ihdms=0; ihdms<ndms; ++ihdms) {
    h_DMS[ihdms] = new TH2F(Form("h_dms_%d",ihdms),h_DMSLabels[ihdms].Data(),nbins_s1_dms,s1_min_dms,s1_max_dms,nbins_f90_dms,f90_min_dms,f90_max_dms);
    h_DMS[ihdms]->GetXaxis()->SetTitle("S1 [PE]");
    h_DMS[ihdms]->GetYaxis()->SetTitle("f_{90}");
  }

 TH2F* h_NR_S1_f90 = new TH2F("h_NR_S1_f90","Nuclear recoil events from AmBe", 200, 0., 1000., 100., 0., 1.);
 h_NR_S1_f90 -> GetXaxis() -> SetTitle("S1 [PE]");
 h_NR_S1_f90 -> GetYaxis() -> SetTitle("f90");


  // S1 vs. tdrift
  outfile->cd(type);
  TDirectory* s1tdrift_dir = gDirectory->mkdir("s1tdrift");
  s1tdrift_dir->cd();
  const Int_t n_s1_tdrift = 2;
  TH2F* h_s1_tdrift[n_s1_tdrift];
  const TString s_s1_tdrift[n_s1_tdrift] = { "No veto cuts", "All cuts"};
  for (Int_t i=0; i<n_s1_tdrift; ++i) {
    h_s1_tdrift[i] = new TH2F(Form("h_s1_tdrift_%d", i), s_s1_tdrift[i], 2000, 0., 10.E+3, 400, 0., 400.); //5 PE S1 binning, 1 us tdrift binning
    h_s1_tdrift[i]->GetXaxis()->SetTitle("S1 [PE]");
    h_s1_tdrift[i]->GetYaxis()->SetTitle("t_{drift} [#mus]");
  }

  // S2/S1 vs. S1
  outfile->cd(type);
  TDirectory* s1_s2s1_dir = gDirectory->mkdir("s1_s2s1");
  s1_s2s1_dir->cd();
  const Int_t n_s1_s2s1 = 1;
  TH2F* h_s1_s2s1[n_s1_s2s1];
  const TString s_s1_s2s1[n_s1_s2s1] = {"No veto cuts, deep R"};
  for (Int_t i=0; i<n_s1_s2s1; ++i) {
    h_s1_s2s1[i] = new TH2F(Form("h_s1_s2s1_%d", i), s_s1_s2s1[i], 200, 0, 1000, 250, 0, 2.5);
    h_s1_s2s1[i]->GetXaxis()->SetTitle("S1 [PE]");
    h_s1_s2s1[i]->GetYaxis()->SetTitle("Log_{10}(S2/S1)");
  }
  
  
  // S2/S1 vs. F90
  outfile->cd(type);
  TDirectory* f90_s2s1_dir = gDirectory->mkdir("f90_s2s1");
  f90_s2s1_dir->cd();
  const Int_t n_f90_s2s1 = 3;
  TH2F* h_f90_s2s1[n_f90_s2s1];
  const TString s_f90_s2s1[n_f90_s2s1] = {"No veto cuts, deep R", "All cuts", "All cuts, deep R"};
  for (Int_t i=0; i<n_f90_s2s1; ++i) {
    h_f90_s2s1[i] = new TH2F(Form("h_f90_s2s1_%d", i), s_f90_s2s1[i], 500, 0, 1, 250, 0, 2.5);
    h_f90_s2s1[i]->GetXaxis()->SetTitle("f_{90}");
    h_f90_s2s1[i]->GetYaxis()->SetTitle("Log_{10}(S2/S1)");
  }

  /*//UNCOMMENT THIS BLOCK IF YOU WANT TO LOOK AT INTEGRATED PULSES OF SEPARATE EVENTS
  outfile->cd(type);                        // FOR S1 WAVEFORMS OF EVENTS
  TDirectory* s1_waveform = gDirectory->mkdir("s1_waveform");
  s1_waveform->cd();
  const Int_t n_s1_waveform = 14;
  TH2F* h_waveformS1[n_s1_waveform];
  Double_t xrange[169];
  
  for(Int_t h=1;h<170;h++){
      xrange[0]=0;
        if(h<=25) xrange[h]=xrange[h-1] + 0.004;
        if(h>25 && h<=34) xrange[h]=xrange[h-1] + 0.1;
        if(h>34 && h<=54) xrange[h]=xrange[h-1] + 0.2;
        if(h>54 && h<=64) xrange[h]=xrange[h-1] + 0.5;
        if(h>64 && h<=74) xrange[h]=xrange[h-1] + 1.;
        if(h>74 && h<=89) xrange[h]=xrange[h-1] + 2.;
        if(h>89 && h<=169) xrange[h]=xrange[h-1] + 5.;
}

  for(Int_t nhs1 = 0; nhs1 < n_s1_waveform; ++nhs1){ 
      h_waveformS1[nhs1] = new TH2F(Form("s1_waveform_%d",nhs1), "", 167, xrange, 1000, 0., 2000.);
      h_waveformS1[nhs1] -> GetXaxis() -> SetTitle("time [ms]");
      h_waveformS1[nhs1] -> GetYaxis() -> SetTitle("S1 integrated, a. u.");
}
  std::cout<< "Histograms defined" << std::endl;
*/

ofstream s1pulse("s1_pulse.txt");
  
  outfile->cd(type);

  
  //-------------------------//
  //     MAIN EVENT LOOP     //
  //-------------------------//
  Int_t tpc_nevents = tpc_events->GetEntries(); // SHOULD BE UNCOMMENTED FOR ALL EVENTS
//   Int_t event_number[14] = {4010, 4162, 4226, 4242, 4261, 4335, 4390, 4502, 4534, 4552, 3011, 3012, 3013, 3014}; // for separate events only
  //if (aar) tpc_nevents = 1e6;
  cout << "Total events: " << tpc_nevents << '\n';
  e.muon_dt = 999;
  int NR_event_counter = 0;
  Double_t livetime_before_12638 = 0;
  Double_t livetime_after_12638 = 0;
  for (Int_t n = 0; n<tpc_nevents; n++) { // for separate events n -> num, tpc_nevents -> size of event_number array
//       Int_t n = event_number[num]; // should be uncommented for separate events
    if (!(n%1000000)) cout << "Processing event " << n << ", " << Int_t(100.*n/tpc_nevents) << "% completed" << endl;
    tpc_events->GetEntry(n);    
//   std::cout<<"Event_number "<<n<<std::endl;  
    if (e.tpc_run_id < 12638) livetime_before_12638 += e.live_time;
    else livetime_after_12638 += e.live_time;
//	cout << "n event=" << n << endl;
//	cout << "vector size " << e.veto_cluster_dtprompt_vec->size() << endl;
//	cout <<  "e.veto_cluster_charge_vec->at(0)" << e.veto_cluster_charge_vec->at(0) << endl;
//	for (unsigned i=0; i<e.veto_cluster_dtprompt_vec->size(); i++){
//	if (e.veto_cluster_dtprompt_vec->at(i)>-0.05 && e.veto_cluster_dtprompt_vec->at(i)<-0.04)
 // 	cout << e.veto_cluster_dtprompt_vec->at(i) << endl;
//	}

	//if (e.veto_cluster_dtprompt_vec->at(0) > -0.05 && e.veto_cluster_dtprompt_vec->at(0) < -0.04 && e.veto_cluster_charge_vec->at(0) > 2400 && e.veto_cluster_charge_vec->at(0) < 3800){	

    // Calculate some variables on the fly
    e.s1_prompt       = e.total_f90 * e.total_s1;
    e.s1_late         = e.total_s1 - e.s1_prompt;
    e.x               = (e.x_masas > -99. ? e.x_masas : (e.x_jasons > -99. ? e.x_jasons : -998.)); //Use Masa's xy by default. If not good, try Jason's. If still not good, junk.
    e.y               = (e.y_masas > -99. ? e.y_masas : (e.y_jasons > -99. ? e.y_jasons : -998.));
    e.r               = (e.r_masas > -99. ? e.r_masas : (e.r_jasons > -99. ? e.r_jasons : -998.));
    e.xycorr_factor   = (e.xycorr_factor_masas > -99. ? e.xycorr_factor_masas : (e.xycorr_factor_jasons > -99. ? e.xycorr_factor_jasons : -998.));
    e.total_s2_xycorr = e.total_s2 * e.xycorr_factor;
//    Int_t s1mf_xbin   = h_s1mf_thresholds->GetXaxis()->FindBin(e.total_s1_corr);
    Int_t s1mf_xbin   = h_s1mf_thresholds->GetXaxis()->FindBin(e.s1_prompt); // S1pmf cut uses s1_prompt binning
    Int_t s1mf_ybin   = h_s1mf_thresholds->GetYaxis()->FindBin(e.tdrift);
    e.s1mf_threshold  = h_s1mf_thresholds->GetBinContent(s1mf_xbin, s1mf_ybin);
    
    // Copied from addvar.C of 50-day analysis
    //e.nr_median_s2    = 5.248/(TMath::Power(0.001102 * e.total_s1_corr,0.7648) + 0.01579) * e.total_s1_corr;
    //e.nr_rms_s2       = -51.06 + 0.2176 * e.nr_median_s2 + 1.321E-5 * TMath::Power(e.nr_median_s2,2.) - 1.150E-9 * TMath::Power(e.nr_median_s2,3.);
    //e.nr_sigma_s2     = (e.total_s2_xycorr - e.nr_median_s2)/e.nr_rms_s2;

    // New S2 over s1
    e.nr_s2overs1 = 4.47255/(TMath::Power(1.27616e-3*e.total_s1_corr, 7.88596e-1) + 1.62597e-2);
    e.s2overs1 = (e.total_s2_xycorr > 0 ? e.total_s2_xycorr/e.total_s1_corr : -999);
    
    e.p0f5000         = (e.npulses > 0 ? e.pulse_us_integrals[4] / e.pulse_us_integrals[6] : -1.);
    
    // Generate cuts.
    // Use booleans so that each cut is evaluated only once per event.
    Bool_t is_prescaled       = isPrescaled(e);
    Bool_t CX_nchan           = cx_nchan(e); // CX#1
    Bool_t CX_baseline        = cx_baseline(e); // CX#2
    Bool_t CX_event_dt        = cx_event_dt(e); // CX#3
    Bool_t CX_file_io         = (aar ? cx_file_io(e) : true); // CX#4 //apply if type is aar. otherwise, don't apply
    Bool_t CX_veto_present    = ((aar || uar) && cx_veto_present(e)) || (ambe);
    Bool_t CX_veto_prompt     = (aar && cx_old_veto_prompt(e))  || (uar && cx_veto_prompt(e))    || (ambe);
    Bool_t CX_veto_delayed    = (aar && cx_old_veto_delayed(e)) || (uar && cx_veto_delayed(e))   || (ambe);
    Bool_t CX_veto_preprompt  = (aar)                           || (uar && cx_veto_preprompt(e)) || (ambe);
    Bool_t CX_veto_muon       = (aar)                           || (uar && cx_veto_muon(e))      || (ambe);
    Bool_t CX_veto_cosmogenic = (aar)                           || (uar && cx_veto_cosmogenic(e))|| (ambe);
    Bool_t CX_single_scatter  = cx_single_scatter(e); // CX#8
    Bool_t CX_trg_time        = cx_trg_time(e); // CX#9
    Bool_t CX_s1_sat          = cx_s1_sat(e); // CX#10
    Bool_t CX_s1_mf           = cx_s1_mf(e); // CX#11
    Bool_t CX_s2_f90          = cx_s2_f90(e); // CX#12
    Bool_t CX_s2_size         = cx_s2_size(e); // CX#13
    Bool_t CX_max_s2          = cx_max_s2(e); // CX#??
    Bool_t CX_s1_range        = cx_s1_range(e); // CX#14
    Bool_t CX_fiducial        = cx_fiducial(e); // CX#15
    Bool_t CX_r16             = cx_r(e, 16.); // CX#16
    Bool_t CX_s2overs1        = cx_s2overs1(e); // CX#17
    Bool_t CX_r10             = cx_r(e, 10.); // CX#18
    Bool_t CX_f90_fast        = cx_f90_fast(e); // CX#19
    Bool_t CX_core            = cx_core(e);
    Bool_t CX_core2           = cx_core2(e);
    Bool_t CX_full_volume     = cx_full_volume(e);

    Bool_t CX_kr85m_1p        = cx_kr85m_1p(e);
    Bool_t CX_kr85m_2p        = cx_kr85m_2p(e);
    Bool_t CX_kr85m_f5000     = cx_kr85m_f5000(e);
    Bool_t CX_kr85m_box       = cx_kr85m_box(e);
    e.kr85mpeaktime = ((CX_kr85m_1p ? e.pulse_peak_time[0] - e.pulse_start_time[0] : 0) +
                       (CX_kr85m_2p ? e.pulse_peak_time[(e.pulse_peak_amp[0] >= e.pulse_peak_amp[1] ? 0 : 1)] - e.pulse_start_time[0] : 0));
    Bool_t CX_kr85m_peaktime  = cx_kr85m_peaktime(e);

    Bool_t CX_f90_s2s1_s1range = cx_f90_s2s1_s1range(e);
    
    // Basics cuts
    Bool_t CX_quality        = (CX_nchan && CX_baseline && CX_event_dt && CX_file_io && CX_veto_present); // CX#QUALITY
    Bool_t CX_veto           = (CX_veto_present && CX_veto_prompt && CX_veto_delayed && CX_veto_preprompt && CX_veto_muon && CX_veto_cosmogenic);
    Bool_t CX_physics        = (CX_quality && CX_veto &&
                                CX_single_scatter && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size &&
                                CX_s1_range && CX_fiducial); // CX#PHYSICS
    Bool_t CX_physics_v1     = (CX_quality && CX_veto &&
                                CX_single_scatter && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size &&
                                CX_fiducial); // CX#PHYSICS-V1
    Bool_t CX_spc            = (CX_physics && CX_r16 && CX_s2overs1); // CX#SPC
    Bool_t CX_test           = (CX_spc && CX_r10 && CX_f90_fast); // CX#TEST
    Bool_t CX_spc_v1         = (CX_physics_v1 && CX_r16 && CX_s2overs1); // CX#SPC-V1
    Bool_t CX_test_v1        = (CX_spc_v1 && CX_r10 && CX_f90_fast); // CX#TEST-V1

    Bool_t isNR              = ambe && (e.total_f90 > ambe_acc99->Eval(e.total_s1_corr)) && (e.veto_roi_lsv_charge_vec->at(0) > 2400.) && (e.veto_roi_lsv_charge_vec->at(0) < 3600.);

    // Quantities
    Float_t s1               = e.total_s1_corr;
    Float_t f90              = e.total_f90;
    Float_t s2               = e.total_s2 * e.xycorr_factor;
    Float_t log10s2overs1    = (s2>0 ? TMath::Log10(s2/s1) : -999);
    
	 
	    // Keep track of time since last muon
	    if (!CX_veto_muon) e.muon_dt = 0.;
    else               e.muon_dt += (e.live_time + e.inhibit_time); // BUG FIX (G.Koh): used to be live + acqui.
	if (e.tdrift > 40. && e.tdrift < 336.){
//	 unsigned i=0;
  //      if (e.veto_cluster_charge_vec->at(i)>2400 && e.veto_cluster_charge_vec->at(i)<3600){

	
	for (unsigned i=0; i<e.veto_cluster_dtprompt_vec->size(); i++){
     if (e.veto_cluster_dtprompt_vec->at(i)>-0.05 && e.veto_cluster_dtprompt_vec->at(i)<-0.04){
     //cout << e.veto_cluster_dtprompt_vec->at(i) << endl;
	 for (unsigned j=0; j<e.veto_cluster_charge_vec->size(); j++){
     if (e.veto_cluster_charge_vec->at(j)>2400 && e.veto_cluster_charge_vec->at(j)<3600){
    
 
    // filling histos
    // h_EventCounter
    h_EventCounter->Fill(NOCUTS);
    if (CX_nchan)                                                             h_EventCounter->Fill(NCHANNEL);
    if (CX_nchan && CX_baseline)                                              h_EventCounter->Fill(BASELINE);
    if (CX_nchan && CX_baseline && CX_event_dt)                               h_EventCounter->Fill(LIVEINHIBIT);
    if (CX_nchan && CX_baseline && CX_event_dt && CX_file_io)                 h_EventCounter->Fill(LONGWAIT);
    if (CX_quality)                                                           h_EventCounter->Fill(VETO_PRESENT);
    if (CX_quality && CX_veto_prompt)                                         h_EventCounter->Fill(VETO_PROMPT);
    if (CX_quality && CX_veto_prompt && CX_veto_delayed)                      h_EventCounter->Fill(VETO_DELAYED);
    if (CX_quality && CX_veto_prompt && CX_veto_delayed && CX_veto_preprompt) h_EventCounter->Fill(VETO_PREPROMPT);
    if (CX_quality && CX_veto_prompt && CX_veto_delayed && CX_veto_preprompt
        && CX_veto_muon)                                                      h_EventCounter->Fill(VETO_MUON);
    if (CX_quality && CX_veto)                                                h_EventCounter->Fill(VETO_COSMOGENIC);
    if (CX_quality && CX_veto && CX_single_scatter)                           h_EventCounter->Fill(NPULSES);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time)            h_EventCounter->Fill(TRIGGERTIME);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat)                                                         h_EventCounter->Fill(S1SATURATION);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf)                                             h_EventCounter->Fill(S1MAXFRAC);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90)                                h_EventCounter->Fill(S2VALID);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size)                  h_EventCounter->Fill(S2MINIMUM);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial)   h_EventCounter->Fill(Z_FIDUCIAL);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial
        && CX_s1_range)                                                       h_EventCounter->Fill(S1RANGE);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial
        && CX_s1_range && CX_r16)                                             h_EventCounter->Fill(RADIUS16);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial
        && CX_s1_range && CX_r16 && CX_s2overs1)                              h_EventCounter->Fill(S2OVERS1);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial
        && CX_s1_range && CX_r16 && CX_s2overs1 && CX_r10)                    h_EventCounter->Fill(RADIUS10);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial
        && CX_s1_range && CX_r16 && CX_s2overs1 && CX_r10 && CX_f90_fast)     h_EventCounter->Fill(FASTf90);

    // Other histograms
    if (CX_quality && CX_veto && CX_single_scatter)                          h_s1_start    ->Fill(e.s1_start_time);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size /*&& CX_fiducial*/
        && CX_s1_range)                                                      h_tdrift      ->Fill(e.tdrift);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size /*&& CX_fiducial*/
        && CX_s1_range && s1<200.)                                            h_tdrift_lowS1->Fill(e.tdrift);

    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        &&!CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial)  h_s1saturation->Fill(s1);
    if (CX_physics)                                                          h_xy_masas ->Fill(e.x_masas, e.y_masas);
    if (CX_physics)                                                          h_xy_jasons->Fill(e.x_jasons, e.y_jasons);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time 
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size /*&& CX_fiducial*/
        && CX_s1_range)                                                      h_rz_masas ->Fill(e.r_masas*e.r_masas/TPC_RMAX, -e.tdrift);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size /*&& CX_fiducial*/
        && CX_s1_range)                                                      h_rz_jasons->Fill(e.r_jasons*e.r_jasons/TPC_RMAX, -e.tdrift);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size /*&& CX_fiducial*/
        && CX_s1_range)                                                      h_rz       ->Fill(e.r*e.r/TPC_RMAX, -e.tdrift);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time 
        && CX_s1_sat && CX_s1_mf && CX_s2_f90/*&& CX_s2_size*/&& CX_fiducial)h_S2xycorr ->Fill(s1, e.total_s2_xycorr);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90/*&& CX_s2_size*/&& CX_fiducial
        && e.r < 0.)                                                        h_s1_badxy ->Fill(s1);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time 
        && CX_s1_sat && CX_s1_mf && CX_s2_f90/*&& CX_s2_size*/&& CX_fiducial
        && e.r < 0.)                                                        h_s2_badxy ->Fill(e.total_s2);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time 
        && CX_s1_sat && CX_s1_mf && CX_s2_f90/*&& CX_s2_size*/&& CX_fiducial
        && e.r < 0.)                                                        h_tdrift_badxy ->Fill(e.tdrift);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time 
        && CX_s1_sat && CX_s1_mf/*&& CX_s2_f90*/&& CX_s2_size && CX_fiducial
        && CX_s1_range)                                                    h_S2f90->Fill(e.total_s2_f90_fixed);
    if (CX_quality && CX_veto && CX_single_scatter && CX_trg_time
        && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size && CX_fiducial
        && CX_s1_range)                                                    h_rdiff->Fill(e.pulse_r_jasons[1] - e.r_jasons);
    
    // tdrift vs S1 histograms
    // use un-corrected S1.
    // turn off s1 saturation because going to high energy and sat is not strong effect.
    // turn off s1mf cut because cut cuts off at high energy
    if (CX_quality /*&& CX_veto*/&& CX_single_scatter
        && CX_trg_time /*&& CX_s1_sat && CX_s1_mf*/ && CX_s2_f90 && CX_s2_size)  h_s1_tdrift[0]->Fill(e.total_s1, e.tdrift, (is_prescaled ? PRESCALE : 1));
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time /*&& CX_s1_sat && CX_s1_mf*/ && CX_s2_f90 && CX_s2_size)  h_s1_tdrift[1]->Fill(e.total_s1, e.tdrift, (is_prescaled ? PRESCALE : 1));

    // Dark matter search (dms)
    if (CX_quality/*&& CX_veto*/&& CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[0]->Fill(s1,f90);
		
    
    ////////////////////////////////////////////S1 Signal ////////////////////////////////////////////////////////
    
/*    Double_t tS1=0; //time for S1 conversion    UNCOMMENT THIS BLOCK IF YOU WANT TO LOOK AT INTEGRATED PULSES OF SEPARATE EVENTS
    for (int h=0; h<170; h++ ){	
	new_v = e.total_s1_long_us_integrals[h];
	t2->Fill();
	h_s1_waveform->Fill(e.total_s1_long_us_integrals[h]);
	h_s2_waveform->Fill(e.total_s2_long_us_integrals[h]);
        
        h_waveformS1[num]->Fill(tS1, e.total_s1_long_us_integrals[h]);
        if(h<=25) tS1+=0.004;
        if(h>25 && h<=34) tS1+=0.1;
        if(h>34 && h<=54) tS1+=0.2;
        if(h>54 && h<=64) tS1+=0.5;
        if(h>64 && h<=74) tS1+=1.;
        if(h>74 && h<=89) tS1+=2.;
        if(h>89 && h<=169) tS1+=5.;
        }
        h_waveformS1[num]->SetMarkerStyle(20);
    */
    
    if(isNR){
    if (s1pulse.is_open()){ //Write all s1 pulses up to 7 usec into text file 
        s1pulse << e.tpc_event_id << " ";
        for(int i = 0; i<54; i++){
            s1pulse << e.total_s1_long_us_integrals[i] << " ";
        }
        s1pulse << " " << endl;
    }
    else{
        cout << "Can't open the txt file" << endl;
    }
    h_NR_S1_f90 -> Fill(s1, f90);
    NR_event_counter++;
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[1]->Fill(s1,f90);
    if (CX_quality && CX_veto_prompt
        && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[2]->Fill(s1,f90);
    if (CX_quality && CX_veto_delayed && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[3]->Fill(s1,f90);
    if (CX_quality && CX_veto_preprompt && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[4]->Fill(s1,f90);
    if (CX_quality && CX_veto_muon && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[5]->Fill(s1,f90);
    if (CX_quality && CX_veto_cosmogenic && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[6]->Fill(s1,f90);
    if (CX_quality /*&& CX_veto*/ /*&& CX_single_scatter*/
        && CX_trg_time && CX_s1_sat/*&& CX_s1_mf*/&& CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[7]->Fill(s1,f90);
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat/*&& CX_s1_mf*/&& CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial)                                   h_DMS[8]->Fill(s1,f90);
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial && CX_s2overs1 && CX_r10)          h_DMS[9]->Fill(s1,f90);
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial/*&& CX_s2overs1*/ && CX_r10)       h_DMS[10]->Fill(s1,f90);
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range && CX_fiducial&& CX_s2overs1*/ && CX_core2)         h_DMS[11]->Fill(s1,f90);

    // S2/S1 vs. S1
    if (CX_quality /*&& CX_veto*/ && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial && CX_r10)                         h_s1_s2s1[0]->Fill(s1,log10s2overs1);

    if (CX_quality /*&& CX_veto*/ && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial && CX_r10 && CX_f90_s2s1_s1range) h_f90_s2s1[0]->Fill(f90,log10s2overs1);
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial && CX_f90_s2s1_s1range)           h_f90_s2s1[1]->Fill(f90,log10s2overs1);
    if (CX_quality && CX_veto && CX_single_scatter
        && CX_trg_time && CX_s1_sat && CX_s1_mf && CX_s2_f90 && CX_s2_size
        /*&& CX_s1_range*/ && CX_fiducial && CX_r10 && CX_f90_s2s1_s1range) h_f90_s2s1[2]->Fill(f90,log10s2overs1);
    
//     cout<<n<<endl;
     }
  }//loop over events
}}}}
  cout << "Livetime before 12638: "<< livetime_before_12638<<endl;
  cout << "Livetime after 12638: "<<livetime_after_12638<<endl;
  cout << "Number of NR events: " << NR_event_counter << endl;
//h_s1_waveform->Scale(norm, "width");
h_s2_waveform->Scale(norm, "width");
t2->Write("", TObject::kOverwrite);  
s1pulse.close(); //Close the txt file with all pulses written.
// Normalize histograms again by mass
  const Float_t full_mass = TMath::Pi() * TPC_RMAX * TPC_RMAX * (FULL_VOL_TDRIFT_MAX - FULL_VOL_TDRIFT_MIN) * VDRIFT / 10. * LAR_DENSITY / 1000.;
  const Float_t core_mass = TMath::Pi() * CORE_R_MAX * CORE_R_MAX * (CORE_TDRIFT_MAX - CORE_TDRIFT_MIN) * VDRIFT / 10. * LAR_DENSITY / 1000.;

}



//------------------------------------------------------------------------------
// Load SLAD files
TChain* load_files(TString mainfile)
{

  // Load and friend the TTrees
  TChain* tpc_events = new TChain("events");
  tpc_events->Add(mainfile);

  if (load_masas_xy) {
    TString xyfile_masas = mainfile;
    xyfile_masas.Remove(xyfile_masas.Length()-5);
    xyfile_masas += "_masas_xy.root";
    
    TChain* xy_masas = new TChain("masas_xy");
    xy_masas->Add(xyfile_masas);
    tpc_events->AddFriend(xy_masas);
    
    TChain* pulse_xy_masas = new TChain("allpulses_xy");
    pulse_xy_masas->Add(xyfile_masas);
    tpc_events->AddFriend(pulse_xy_masas);
  }

  if (load_xylocator_xy) {
    TString xyfile_jasons = mainfile;
    xyfile_jasons.Remove(xyfile_jasons.Length()-5);
    xyfile_jasons += "_xylocator_xy.root";

    TChain* xy_jasons = new TChain("xylocator_xy");
    xy_jasons->Add(xyfile_jasons);
    tpc_events->AddFriend(xy_jasons);
    
    TChain* pulse_xy_jasons = new TChain("allpulses_xyl_xy");
    pulse_xy_jasons->Add(xyfile_jasons);
    tpc_events->AddFriend(pulse_xy_jasons);
  }

  if (load_aww_xy) {
    TString xyfile_aww = mainfile;
    xyfile_aww.Remove(xyfile_aww.Length()-5);
    xyfile_aww += "_aww_xy.root";
    
    TChain* xy_aww = new TChain("aww_xy");
    xy_aww->Add(xyfile_aww);
    tpc_events->AddFriend(xy_aww);
    
    TChain* pulse_xy_aww = new TChain("allpulses_aww_xy");
    pulse_xy_aww->Add(xyfile_aww);
    tpc_events->AddFriend(pulse_xy_aww);
  }

  if (load_allpulses) {
    TString pulsefile = mainfile;
    pulsefile.Remove(pulsefile.Length()-5);
    pulsefile += "_allpulses.root";
    
    TChain* pulse_info = new TChain("pulse_info");
    pulse_info->Add(pulsefile);
    tpc_events->AddFriend(pulse_info);
  }

  if (load_veto) {
    TString vetofile = mainfile;
    vetofile.Remove(vetofile.Length()-5);
    vetofile += "_veto_cluster.root";
  
    TChain* veto_events = new TChain("veto");
    veto_events->Add(vetofile);
    tpc_events->AddFriend(veto_events);
  }

  return tpc_events;
}


//------------------------------------------------------------------------------
// Main method. Load DST files and invoke event_loop().
void data_validation() {
  // Prevent canvases from being drawn.
  gROOT->SetBatch(kTRUE);
  gROOT->Reset();
  SetMyStyle();

  TStopwatch* clock = new TStopwatch();
  clock->Start();

  cout << "Saving output to " << outfile->GetName() << '\n';

  Bool_t doUAr     = 1;
  
  if (doUAr) {
    load_allpulses = true;
    load_masas_xy = true;
    load_xylocator_xy = false;
    load_veto = true; // CHANGED IT FROM TRUE
    
    // The main SLAD file containing the data we want.
//     TString mainfile = "/home/agr/AGR/DarkSide/Run020644_SLAD.root";
    TString mainfile = "/home/agr/AGR/DarkSide/AmBe_160nps_SLAD_v3_3_0_merged.root";
    TChain* tpc_events = load_files(mainfile);
    event_loop(tpc_events, "ambe"); //CHANGED from ambe
  }

  

  outfile->Write();
  
  outfile->Close();
  
  cout << "Done! " << clock->RealTime() << " s." << endl;
}
