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
 */

#ifndef PATTERN_TILING_H_
#define PATTERN_TILING_H_

#include "../utils/commons.h"

/*
 * Pattern Tiling Generator
 */
typedef struct {
  // Current tile dimensions
  uint64_t tile_offset;
  uint64_t tile_wide;
  uint64_t tile_tall;
  uint64_t tile_next_offset_inc;
  // Complete tile dimensions (defaults)
  uint64_t pattern_band_width;
  uint64_t pattern_tile_offset;
  uint64_t pattern_tile_wide;
  uint64_t pattern_tile_tall;
  uint64_t pattern_max_error;
  uint64_t pattern_remaining_length;
  uint64_t sequence_length;
} pattern_tiling_t;

/*
 * Pattern tiling Generation
 */
void pattern_tiling_init(
    pattern_tiling_t* const pattern_tiling,
    const uint64_t pattern_length,
    const uint64_t pattern_tile_length,
    const uint64_t sequence_length,
    const uint64_t max_error);
void pattern_tiling_next(
    pattern_tiling_t* const pattern_tiling);

#endif /* PATTERN_TILING_H_ */
