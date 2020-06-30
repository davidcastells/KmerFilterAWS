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
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Edit-Distance based BPM alignment algorithm
 */

#ifndef EDIT_BPM_H_
#define EDIT_BPM_H_

#include "../utils/commons.h"
#include "edit_column.h"
#include "../system/mm_allocator.h"

/*
 * BPM Pattern
 */
typedef struct {
  /* BMP Pattern */
  char* pattern;                // Raw pattern
  uint64_t* PEQ;                // Pattern equalities (Bit vector for Myers-DP)
  uint64_t pattern_length;      // Length
  uint64_t pattern_num_words64; // ceil(Length / |w|)
  uint64_t pattern_mod;         // Length % |w|
  /* BPM Auxiliary data */
  uint64_t* P;
  uint64_t* M;
  uint64_t* level_mask;
  int64_t* score;
  int64_t* init_score;
  uint64_t* pattern_left;
} bpm_pattern_t;

/*
 * Setup
 */
void edit_bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,
    char* const pattern,
    const int pattern_length,
    mm_allocator_t* const mm_allocator);
void edit_bpm_pattern_free(
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator);

/*
 * Edit distance computation using BPM
 */
int edit_bpm_distance_compute(
    bpm_pattern_t* const bpm_pattern,
    char* const text,
    const int text_length);
int edit_bpm_distance_compute_cutoff(
    bpm_pattern_t* const bpm_pattern,
    char* const text,
    const int text_length,
    int max_distance,
    const bool quick_abandon);

#endif /* EDIT_BPM_H_ */
