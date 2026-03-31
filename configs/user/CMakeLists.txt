cmake_minimum_required(VERSION 3.20 FATAL_ERROR) #minimum version of CMAKE required
project(UserActions)

# Load ActRoot custom script
set(ACTROOT $ENV{ACTROOT})
include(${ACTROOT}/Scripts.cmake)

#And compile
add_userlibrary(NAME RecRANSAC SOURCES RecRANSAC.h RecRANSAC.cxx LINK ActAlgorithm)
add_userlibrary(NAME FilterDecay SOURCES FilterDecay.h FilterDecay.cxx LINK ActAlgorithm)
add_userlibrary(NAME ConfMergerPerRun SOURCES ConfMergerPerRun.h ConfMergerPerRun.cxx LINK ActAlgorithm)
add_userlibrary(NAME PCASplit SOURCES PCASplit.h PCASplit.cxx LINK ActAlgorithm)
