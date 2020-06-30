/*
 *  Wavefront Alignments Algorithms
 *  Copyright (c) 2020 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of Wavefront Alignments Algorithms.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: Fast Mapping-Candidates Filtering Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Fast Mapping-Candidates Filtering Benchmark
 */

#include "../utils/commons.h"
#include "../utils/input_text.h"
#include "../system/profiler_timer.h"
#include "../benchmark/benchmark_utils.h"
#include "../benchmark/benchmark_edit_alg.h"
#include "../benchmark/benchmark_kmer_filter.h"
#include "FPGAKmerFilter.h"

/*
 * Algorithms
 */
typedef enum {
  filter_edit_dp,
  filter_edit_bpm,
  filter_kmer_nway,
  filter_kmer_fpga,
} filter_type;

/*
 * Parameters
 */
typedef struct {
  // Input
  char *algorithm;
  char *input;
  float max_error;
  // Specifics
  float bandwidth;
  int kmer_length;
  // Profile
  profiler_timer_t timer_global;
  int progress;
  // Misc
  bool check;
  bool verbose;
} benchmark_args;

benchmark_args parameters;



int init_parameters()
{
  // Input
  parameters.algorithm=NULL;
  parameters.input=NULL;
  parameters.max_error=0.05;
  // Specifics
  parameters.bandwidth = -1.0;
  parameters.kmer_length = 5;
  // Profile
  parameters.progress = 100000;
  // Misc
  parameters.check = false;
  parameters.verbose = false;
          
  return 0;
}

int dummy = init_parameters();

/*
 * UTest
 */
void filter_test() {
//  // Patters & Texts
//  char* pattern = "ACGTACGT";
//  char* text    = "TGCATGCA";
//
//  // Configure filter_input
//  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
//
//  // Free
//  mm_allocator_delete(mm_allocator);
}
/*
 * Benchmark
 */
void filter_benchmark(const filter_type filter) 
{
  // Parameters
  FILE *input_file = NULL;
  char *line1 = NULL, *line2 = NULL;
  int line1_length=0, line2_length=0;
  size_t line1_allocated=0, line2_allocated=0;
  filter_input_t filter_input;
  // Init
  timer_restart(&(parameters.timer_global));
  input_file = fopen(parameters.input, "r");
  
  if (input_file==NULL)
  {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input);
    exit(1);
  }
  filter_input_clear(&filter_input);
  filter_input.check = parameters.check;
  filter_input.verbose = parameters.verbose;
  filter_input.mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  // Read-filter loop
  int seq_processed = 0, progress = 0;
  
  FPGAKmerFilter fpga;
        fpga.initOpenCL(0);
      fpga.initKernels(1, "fpga" );

  timer_reset(&filter_input.timer);

  
  while (true) 
  {
    // Read queries
    line1_length = getline(&line1,&line1_allocated,input_file);
    if (line1_length==-1) break;
    line2_length = getline(&line2,&line2_allocated,input_file);
    if (line1_length==-1) break;
    // Configure input
    filter_input.sequence_id = seq_processed;
    filter_input.pattern = line1+1;
    filter_input.pattern_length = line1_length-2;
    filter_input.pattern[filter_input.pattern_length] = '\0';
    filter_input.text = line2+1;
    filter_input.text_length = line2_length-2;
    filter_input.text[filter_input.text_length] = '\0';
    if (parameters.max_error >= 1.0) {
      filter_input.max_error = parameters.max_error;
    } else {
      filter_input.max_error = ceil(parameters.max_error*filter_input.pattern_length);
    }
    int bandwidth;
    if (parameters.bandwidth >= 1.0) {
      bandwidth = parameters.bandwidth;
    } else if (parameters.bandwidth >= 0.0) {
      bandwidth = ceil(parameters.bandwidth*filter_input.pattern_length);
    } else {
      bandwidth = -1;
    }
    // Filter
    switch (filter) {
      case filter_edit_dp:
        benchmark_edit_dp(&filter_input,bandwidth);
        break;
      case filter_edit_bpm:
        benchmark_edit_bpm(&filter_input,bandwidth);
        break;
      case filter_kmer_nway:
        benchmark_kmer_filter(&filter_input,parameters.kmer_length);
        break;
      case filter_kmer_fpga:
        fpga.addInput(&filter_input,parameters.kmer_length);
        break;
      default:
        fprintf(stderr,"Algorithm unknown or not implemented\n");
        exit(1);
        break;
    }
    // Update progress
    ++seq_processed;
    // DEBUG mm_allocator_print(stderr,filter_input.mm_allocator,true);
    if (++progress == parameters.progress) 
    {
      progress = 0;
      // Compute statistics
      const uint64_t time_elapsed = timer_elapsed_ns(&(filter_input.timer));
      const float time_filter_rate = (float)seq_processed/(float)TIMER_CONVERT_NS_TO_S(time_elapsed);
      fprintf(stderr,"...processed %d sequences (filter=%2.3f sequences/s)\n",seq_processed,time_filter_rate);
    }
  }
  
  if (filter == filter_kmer_fpga)
  {
      //fpga.initKernels(1, "emulator");
	fpga.m_verbose = true;
      fpga.computeAll();
      //fpga.destroy();
      
      fpga.finalizeKernels();
      fpga.finalizeOpenCL();
  }
  
  timer_stop(&(parameters.timer_global));
  // Print benchmark results
  fprintf(stderr,"[Benchmark]\n");
  fprintf(stderr,"=> Total.sequences        %d\n",seq_processed);
  fprintf(stderr,"=> Time.Benchmark      ");
  timer_print(stderr,&parameters.timer_global,NULL);
  fprintf(stderr,"  => Time.Filter       ");
  timer_print(stderr,&filter_input.timer,&parameters.timer_global);
  if (parameters.check) {
    fprintf(stderr,"=> Check\n");
    fprintf(stderr,"  => TP.Hit       %d (%2.3f)\n",(filter_input.candidates_tp+filter_input.candidates_tn),
        100.0f*(float)(filter_input.candidates_tp+filter_input.candidates_tn)/(float)filter_input.candidates_total);
    fprintf(stderr,"    => TP.OK.IN   %d (%2.3f)\n",filter_input.candidates_tp,
        100.0f*(float)filter_input.candidates_tp/(float)filter_input.candidates_total);
    fprintf(stderr,"    => TN.OK.OUT  %d (%2.3f)\n",filter_input.candidates_tn,
        100.0f*(float)filter_input.candidates_tn/(float)filter_input.candidates_total);
    fprintf(stderr,"  => FP.Noise     %d (%2.3f)\n",filter_input.candidates_fp,
        100.0f*(float)filter_input.candidates_fp/(float)filter_input.candidates_total);
    fprintf(stderr,"  => FN.Miss      %d (%2.3f)\n",filter_input.candidates_fn,
        100.0f*(float)filter_input.candidates_fn/(float)filter_input.candidates_total);
  }
  // Free
  fclose(input_file);
  mm_allocator_delete(filter_input.mm_allocator);
  free(line1);
  free(line2);
}
/*
 * Generic Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./filters_benchmark -a algorithm -i input                       \n"
      "      Options::                                                      \n"
      "        [Input]                                                      \n"
      "          --algorithm|-a <ALGORITHM>                                 \n"
      "            [edit-filters]                                           \n"
      "              edit-dp                                                \n"
      "              edit-bpm                                               \n"
      "            [kmer-filters]                                           \n"
      "              kmer-filter                                            \n"
      "          --input|-i <FILE>                                          \n"
      "          --max-error|-e <INT>|<FLOAT>       (default=0.05)          \n"
      "        [Specifics]                                                  \n"
      "          --bandwidth|-b <INT>|<FLOAT>       (default=disabled)      \n"
      "          --kmer-length|-k [3..7]            (default=5)             \n"
      "        [Misc]                                                       \n"
      "          --progress|-P <INT>                                        \n"
      "          --help|-h                                                  \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* Input */
    { "algorithm", required_argument, 0, 'a' },
    { "input", required_argument, 0, 'i' },
    { "max-error", required_argument, 0, 'e' },
    /* Specifics */
    { "bandwidth", required_argument, 0, 'b' },
    { "kmer-length", required_argument, 0, 'k' },
    /* Misc */
    { "progress", required_argument, 0, 'P' },
    { "check", no_argument, 0, 'c' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"a:i:e:b:k:P:cvh",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * Input
     */
    case 'a':
      parameters.algorithm = optarg;
      break;
    case 'i':
      parameters.input = optarg;
      break;
    case 'e':
      parameters.max_error = atof(optarg);
      break;
    /*
     * Specific parameters
     */
    case 'b': // --bandwidth
      parameters.bandwidth = atof(optarg);
      break;
    case 'k': // --kmer-length
      parameters.kmer_length = atoi(optarg);
      break;
    /*
     * Misc
     */
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'c':
      parameters.check = true;
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
}
int main(int argc,char* argv[]) {
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Select option
  if (strcmp(parameters.algorithm,"test")==0) {
    filter_test();
  /* Edit */
  } else if (strcmp(parameters.algorithm,"edit-dp")==0) {
    filter_benchmark(filter_edit_dp);
  } else if (strcmp(parameters.algorithm,"edit-bpm")==0) {
    filter_benchmark(filter_edit_bpm);
  } else if (strcmp(parameters.algorithm,"kmer-filter")==0) {
    filter_benchmark(filter_kmer_nway);
  } else if (strcmp(parameters.algorithm,"kmer-fpga")==0) {
    filter_benchmark(filter_kmer_fpga);
  } else {
    fprintf(stderr,"Algorithm '%s' not recognized\n",parameters.algorithm);
    exit(1);
  }
}



