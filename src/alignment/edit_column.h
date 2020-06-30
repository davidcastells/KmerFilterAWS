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
 * DESCRIPTION: Dynamic-programming alignment table
 */

#ifndef EDIT_COLUMN_H_
#define EDIT_COLUMN_H_

#include "../utils/commons.h"
#include "../system/mm_allocator.h"

/*
 * Constants
 */
#define SCORE_MAX (10000000)

/*
 * DP Table
 */
typedef struct {
  // Dimensions
  int num_rows;
  int num_columns;
  // DP Columns (current and previous)
  int* column_curr;
  int* column_prev;
} edit_column_t;

/*
 * Setup
 */
void edit_column_allocate(
    edit_column_t* const edit_column,
    const int pattern_length,
    const int text_length,
    mm_allocator_t* const mm_allocator);
void edit_column_free(
    edit_column_t* const edit_column,
    mm_allocator_t* const mm_allocator);

#endif /* EDIT_COLUMN_H_ */
