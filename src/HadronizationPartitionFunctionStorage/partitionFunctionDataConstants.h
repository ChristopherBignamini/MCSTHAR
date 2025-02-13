#ifndef PARTITIONFUNCTIONDATACONSTANTS_H
#define PARTITIONFUNCTIONDATACONSTANTS_H

#include <string>

using namespace std;

/**
* @brief Partition function data/archive/summary constant file names (partial) strings
*
* @author Christopher Bignamini
*/

/**
* Partition function archive file (relative) name
*/
// TODO: static or not?
// TODO: unify s_ everywhere
static const string s_partitionFunctionArchiveFileName("partition_function_archive_data.xml");

/**
* Partition function data summary file (relative) name
*/
static const string s_partitionFunctionDataSummaryFileName("partition_function_data_summary.xml");

/**
* Partition function data file name prefix
*/
static const string s_partitionFunctionDataFileNamePrefix("partition_function_data_");

#endif