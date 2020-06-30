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
 * DESCRIPTION: Dynamic-programming alignment algorithm for Levenshtein distance (edit)
 */

#ifndef EDIT_DP_H_
#define EDIT_DP_H_

#include "../utils/commons.h"
#include "../alignment/edit_column.h"

/*
 * Edit distance computation [Ends-free] [DP-Columns]
 */
int edit_dp_distance(
    edit_column_t* const edit_column,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length);
/*
 * Edit distance computation [Ends-free] [DP-Columns] [Banded]
 */
int edit_dp_distance_banded(
    edit_column_t* const edit_column,
    const char* const pattern,
    const int pattern_length,
    const char* const text,
    const int text_length,
    const int bandwidth);

#endif /* EDIT_DP_H_ */
