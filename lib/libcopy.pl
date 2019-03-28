#!/usr/bin/perl

$dir_target = "/gpfs/fs03/t2k/beam/work/nchikuma/ROOT/v5r24p00n02/amd64_linux26/lib/root";
$dir_source = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/lib";

if(-d $dir_target and -d $dir_source){
  print "Target directory = $dir_target\n";
  print "Source directory = $dir_source\n";

  system "ln -sf $dir_source/EVENTSUMMARY.so         $dir_target/libEVENTSUMMARY.so";
  system "ln -sf $dir_source/HitSummary.so           $dir_target/libHitSummary.so";
  system "ln -sf $dir_source/SimHitSummary.so        $dir_target/libSimHitSummary.so";
  system "ln -sf $dir_source/SimVertexSummary.so     $dir_target/libSimVertexSummary.so";
  system "ln -sf $dir_source/SimParticleSummary.so   $dir_target/libSimParticleSummary.so";
  system "ln -sf $dir_source/BeamInfoSummary.so      $dir_target/libBeamInfoSummary.so";
  system "ln -sf $dir_source/BasicReconSummary.so    $dir_target/libBasicReconSummary.so";
  system "ln -sf $dir_source/FirstReducSummary.so    $dir_target/libFirstReducSummary.so";
  system "ln -sf $dir_source/NeutInfoSummary.so      $dir_target/libNeutInfoSummary.so";
  system "ln -sf $dir_source/TrackSummary.so         $dir_target/libTrackSummary.so";
  system "ln -sf $dir_source/TwoDimReconSummary.so   $dir_target/libTwoDimReconSummary.so";
  system "ln -sf $dir_source/ThreeDimReconSummary.so $dir_target/libThreeDimReconSummary.so";
}

#$dir_target = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/INGRID/INGRID/v1r1/src";
#$dir_source = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/lib";
#
#if(-d $dir_target and -d $dir_source){
#
#  print "Target directory = $dir_target\n";
#  print "Source directory = $dir_source\n";
#
#  system "ln -sf $dir_source/EVENTSUMMARY.h         $dir_target/EVENTSUMMARY.h";
#  system "ln -sf $dir_source/HitSummary.h           $dir_target/HitSummary.h";
#  system "ln -sf $dir_source/SimHitSummary.h        $dir_target/SimHitSummary.h";
#  system "ln -sf $dir_source/SimVertexSummary.h     $dir_target/SimVertexSummary.h";
#  system "ln -sf $dir_source/SimParticleSummary.h   $dir_target/SimParticleSummary.h";
#  system "ln -sf $dir_source/BeamInfoSummary.h      $dir_target/BeamInfoSummary.h";
#  system "ln -sf $dir_source/BasicReconSummary.h    $dir_target/BasicReconSummary.h";
#  system "ln -sf $dir_source/FirstReducSummary.h    $dir_target/FirstReducSummary.h";
#  system "ln -sf $dir_source/NeutInfoSummary.h      $dir_target/NeutInfoSummary.h";
#  system "ln -sf $dir_source/TrackSummary.h         $dir_target/TrackSummary.h";
#  system "ln -sf $dir_source/TwoDimReconSummary.h   $dir_target/TwoDimReconSummary.h";
#  system "ln -sf $dir_source/ThreeDimReconSummary.h $dir_target/ThreeDimReconSummary.h";
##
#  system "ln -sf $dir_source/EVENTSUMMARY.cc         $dir_target/EVENTSUMMARY.cc";
#  system "ln -sf $dir_source/HitSummary.cc           $dir_target/HitSummary.cc";
#  system "ln -sf $dir_source/SimHitSummary.cc        $dir_target/SimHitSummary.cc";
#  system "ln -sf $dir_source/SimVertexSummary.cc     $dir_target/SimVertexSummary.cc";
#  system "ln -sf $dir_source/SimParticleSummary.cc   $dir_target/SimParticleSummary.cc";
#  system "ln -sf $dir_source/BeamInfoSummary.cc      $dir_target/BeamInfoSummary.cc";
#  system "ln -sf $dir_source/BasicReconSummary.cc    $dir_target/BasicReconSummary.cc";
#  system "ln -sf $dir_source/FirstReducSummary.cc    $dir_target/FirstReducSummary.cc";
#  system "ln -sf $dir_source/NeutInfoSummary.cc      $dir_target/NeutInfoSummary.cc";
#  system "ln -sf $dir_source/TrackSummary.cc         $dir_target/TrackSummary.cc";
#  system "ln -sf $dir_source/TwoDimReconSummary.cc   $dir_target/TwoDimReconSummary.cc";
#  system "ln -sf $dir_source/ThreeDimReconSummary.cc $dir_target/ThreeDimReconSummary.cc";
#}
