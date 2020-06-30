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

#include "edit_bpm_distance.h"
#include "../utils/dna_text.h"

/*
 * Constants
 */
#define BPM_W64_LENGTH       UINT64_LENGTH
#define BPM_W64_SIZE         UINT64_SIZE
#define BPM_W64_ONES         UINT64_MAX
#define BPM_W64_MASK         (1ull<<63)
#define BPM_ALPHABET_LENGTH  (5)

/*
 * Pattern Accessors
 */
#define BPM_PATTERN_PEQ_IDX(word_pos,encoded_character)   ((word_pos*BPM_ALPHABET_LENGTH)+(encoded_character))
#define BPM_PATTERN_BDP_IDX(position,num_words,word_pos)  ((position)*(num_words)+(word_pos))

/*
 * Advance block functions (Improved)
 *   const @vector Eq,mask;
 *   return (Pv,Mv,PHout,MHout);
 */
#define BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout) \
  /* Computes modulator vector {Xv,Xh} ( cases A&C ) */ \
  const uint64_t Xv = Eq | Mv; \
  const uint64_t _Eq = Eq | MHin; \
  const uint64_t Xh = (((_Eq & Pv) + Pv) ^ Pv) | _Eq; \
  /* Calculate Hout */ \
  uint64_t Ph = Mv | ~(Xh | Pv); \
  uint64_t Mh = Pv & Xh; \
  /* Account Hout that propagates for the next block */ \
  PHout = (Ph & mask)!=0; \
  MHout = (Mh & mask)!=0; \
  /* Hout become the Hin of the next cell */ \
  Ph <<= 1; \
  Mh <<= 1; \
  /* Account Hin coming from the previous block */ \
  Ph |= PHin; \
  Mh |= MHin; \
  /* Finally, generate the Vout */ \
  Pv = Mh | ~(Xv | Ph); \
  Mv = Ph & Xv
/*
 * Setup
 */
void edit_bpm_pattern_compile(
    bpm_pattern_t* const bpm_pattern,
    char* const pattern,
    const int pattern_length,
    mm_allocator_t* const mm_allocator) {
  // Calculate dimensions
  const uint64_t pattern_num_words64 = DIV_CEIL(pattern_length,BPM_W64_LENGTH);
  const uint64_t PEQ_length = pattern_num_words64*BPM_W64_LENGTH;
  const uint64_t pattern_mod = pattern_length%BPM_W64_LENGTH;
  // Init fields
  bpm_pattern->pattern = pattern;
  bpm_pattern->pattern_length = pattern_length;
  bpm_pattern->pattern_num_words64 = pattern_num_words64;
  bpm_pattern->pattern_mod = pattern_mod;
  // Allocate memory
  const uint64_t aux_vector_size = pattern_num_words64*BPM_W64_SIZE;
  const uint64_t PEQ_size = BPM_ALPHABET_LENGTH*aux_vector_size;
  const uint64_t score_size = pattern_num_words64*UINT64_SIZE;
  const uint64_t total_memory = PEQ_size + 3*aux_vector_size + 2*score_size + (pattern_num_words64+1)*UINT64_SIZE;
  void* memory = mm_allocator_malloc(mm_allocator,total_memory);
  bpm_pattern->PEQ = (uint64_t*) memory; 
  memory += PEQ_size;
  bpm_pattern->P = (uint64_t*)memory; 
  memory += aux_vector_size;
  bpm_pattern->M = (uint64_t*)memory; 
  memory += aux_vector_size;
  bpm_pattern->level_mask = (uint64_t*)memory;
  memory += aux_vector_size;
  bpm_pattern->score = (int64_t*)memory; 
  memory += score_size;
  bpm_pattern->init_score =(int64_t*) memory;
  memory += score_size;
  bpm_pattern->pattern_left = (uint64_t*) memory;
  // Init PEQ
  memset(bpm_pattern->PEQ,0,PEQ_size);
  uint64_t i;
  for (i=0;i<pattern_length;++i) {
    const uint8_t enc_char = dna_encode(pattern[i]);
    if (enc_char==ENC_DNA_CHAR_N) continue; // N's Inequality
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,enc_char)] |= mask;
  }
  for (;i<PEQ_length;++i) { // Padding
    const uint64_t block = i/BPM_W64_LENGTH;
    const uint64_t mask = 1ull<<(i%BPM_W64_LENGTH);
    uint64_t j;
    for (j=0;j<BPM_ALPHABET_LENGTH;++j) {
      bpm_pattern->PEQ[BPM_PATTERN_PEQ_IDX(block,j)] |= mask;
    }
  }
  // Init auxiliary data
  uint64_t pattern_left = pattern_length;
  const uint64_t top = pattern_num_words64-1;
  memset(bpm_pattern->level_mask,0,aux_vector_size);
  for (i=0;i<top;++i) {
    bpm_pattern->level_mask[i] = BPM_W64_MASK;
    bpm_pattern->init_score[i] = BPM_W64_LENGTH;
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  for (;i<=pattern_num_words64;++i) {
    bpm_pattern->pattern_left[i] = pattern_left;
    pattern_left = (pattern_left > BPM_W64_LENGTH) ? pattern_left-BPM_W64_LENGTH : 0;
  }
  if (pattern_mod > 0) {
    const uint64_t mask_shift = pattern_mod-1;
    bpm_pattern->level_mask[top] = 1ull<<(mask_shift);
    bpm_pattern->init_score[top] = pattern_mod;
  } else {
    bpm_pattern->level_mask[top] = BPM_W64_MASK;
    bpm_pattern->init_score[top] = BPM_W64_LENGTH;
  }
}
void edit_bpm_pattern_free(
    bpm_pattern_t* const bpm_pattern,
    mm_allocator_t* const mm_allocator) {
  mm_allocator_free(mm_allocator,bpm_pattern->PEQ);
}
/*
 * BPM Auxiliary functions
 */
void bpm_reset_search(
    const uint64_t num_words,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score) {
  // Reset score,P,M
  uint64_t i;
  P[0] = BPM_W64_ONES;
  M[0] = 0;
  score[0] = init_score[0];
  for (i=1;i<num_words;++i) {
    P[i] = BPM_W64_ONES;
    M[i] = 0;
    score[i] = score[i-1] + init_score[i];
  }
}
void bpm_reset_search_cutoff(
    uint8_t* const top_level,
    uint64_t* const P,
    uint64_t* const M,
    int64_t* const score,
    const int64_t* const init_score,
    const uint64_t max_distance) {
  // Calculate the top level (maximum bit-word for cut-off purposes)
  const uint8_t y = (max_distance>0) ? (max_distance+(BPM_W64_LENGTH-1))/BPM_W64_LENGTH : 1;
  *top_level = y;
  // Reset score,P,M
  uint64_t i;
  P[0] = BPM_W64_ONES;
  M[0] = 0;
  score[0] = init_score[0];
  for (i=1;i<y;++i) {
    P[i] = BPM_W64_ONES;
    M[i] = 0;
    score[i] = score[i-1] + init_score[i];
  }
}
/*
 * Advance block functions
 */
int8_t T_hout_64[2][2] = {{0,-1},{1,1}};
uint8_t P_hin_64[3] = {0, 0, 1L};
uint8_t N_hin_64[3] = {1L, 0, 0};
int8_t bpm_advance_block(
    uint64_t Eq,
    const uint64_t mask,
    uint64_t Pv,
    uint64_t Mv,
    const int8_t hin,
    uint64_t* const Pv_out,
    uint64_t* const Mv_out) {
  uint64_t Ph, Mh;
  uint64_t Xv, Xh;
  int8_t hout=0;

  Xv = Eq | Mv;
  Eq |= N_hin_64[hin];
  Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

  Ph = Mv | ~(Xh | Pv);
  Mh = Pv & Xh;

  hout += T_hout_64[(Ph & mask)!=0][(Mh & mask)!=0];

  Ph <<= 1;
  Mh <<= 1;

  Mh |= N_hin_64[hin];
  Ph |= P_hin_64[hin];

  Pv = Mh | ~(Xv | Ph);
  Mv = Ph & Xv;

  *Pv_out=Pv;
  *Mv_out=Mv;

  return hout;
}
/*
 * BPM Distance Compute
 */
int edit_bpm_distance_compute(
    bpm_pattern_t* const bpm_pattern,
    char* const text,
    const int text_length) {
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  // Initialize search
  int min_score = INT_MAX;
  bpm_reset_search(num_words64,P,M,score,init_score);
  // Advance in DP-bit_encoded matrix
  uint64_t text_position;
  for (text_position=0;text_position<text_length;++text_position) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    int8_t carry;
    uint64_t i;
    for (i=0,carry=0;i<num_words64;++i) {
      uint64_t* const Py = P+i;
      uint64_t* const My = M+i;
      carry = bpm_advance_block(PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)],level_mask[i],*Py,*My,carry+1,Py,My);
      score[i] += carry;
    }
    // Check match
    if (score[num_words64-1] < min_score) {
      min_score = score[num_words64-1];
    }
  }
  // Return results
  return min_score;
}
int edit_bpm_distance_compute_cutoff(
    bpm_pattern_t* const bpm_pattern,
    char* const text,
    const int text_length,
    int max_distance,
    const bool quick_abandon) {
  // Pattern variables
  const uint64_t* PEQ = bpm_pattern->PEQ;
  const uint64_t num_words64 = bpm_pattern->pattern_num_words64;
  uint64_t* const P = bpm_pattern->P;
  uint64_t* const M = bpm_pattern->M;
  const uint64_t* const level_mask = bpm_pattern->level_mask;
  int64_t* const score = bpm_pattern->score;
  const int64_t* const init_score = bpm_pattern->init_score;
  const uint64_t* const pattern_left = bpm_pattern->pattern_left;
  // Initialize search
  if (max_distance >= bpm_pattern->pattern_length) {
    max_distance = bpm_pattern->pattern_length-1; // Correct max-distance
  }
  const uint64_t max_distance__1 = max_distance+1;
  const uint8_t top = num_words64-1;
  uint8_t top_level;
  int min_score = INT_MAX;
  bpm_reset_search_cutoff(&top_level,P,M,score,init_score,max_distance);
  // Advance in DP-bit_encoded matrix
  uint64_t text_position, text_left=text_length;
  for (text_position=0;text_position<text_length;++text_position,--text_left) {
    // Fetch next character
    const uint8_t enc_char = dna_encode(text[text_position]);
    // Advance all blocks
    uint64_t i,PHin=0,MHin=0,PHout,MHout;
    for (i=0;i<top_level;++i) {
      uint64_t Pv = P[i];
      uint64_t Mv = M[i];
      const uint64_t mask = level_mask[i];
      const uint64_t Eq = PEQ[BPM_PATTERN_PEQ_IDX(i,enc_char)];
      /* Compute Block */
      BPM_ADVANCE_BLOCK(Eq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
      /* Save Block Pv,Mv */
      P[i]=Pv;
      M[i]=Mv;
      /* Adjust score and swap propagate Hv */
      score[i] += PHout-MHout;
      PHin=PHout;
      MHin=MHout;
    }
    // Cut-off
    const uint8_t last = top_level-1;
    if (score[last]<=max_distance__1 && last<top) {
      const uint64_t last_score = score[last]+(MHin-PHin);
      const uint64_t Peq = PEQ[BPM_PATTERN_PEQ_IDX(top_level,enc_char)];
      if (last_score<=max_distance && (MHin || (Peq & 1))) {
        // Init block V
        uint64_t Pv = BPM_W64_ONES;
        uint64_t Mv = 0;
        const uint64_t mask = level_mask[top_level];
        /* Compute Block */
        BPM_ADVANCE_BLOCK(Peq,mask,Pv,Mv,PHin,MHin,PHout,MHout);
        /* Save Block Pv,Mv */
        P[top_level]=Pv;
        M[top_level]=Mv;
        /* Set score & increment the top level block */
        score[top_level] = last_score + init_score[top_level] + (PHout-MHout);
        ++top_level;
      } else {
        while (score[top_level-1] > (max_distance+init_score[top_level-1])) {
          --top_level;
        }
      }
    } else {
      while (score[top_level-1] > (max_distance+init_score[top_level-1])) {
        --top_level;
      }
    }
    // Check match
    const int64_t current_score = score[top_level-1];
    if (top_level==num_words64 && current_score<=max_distance) {
      if (current_score < min_score) min_score = current_score;
    } else if (quick_abandon) {
      // Quick abandon, if it doesn't match (bounded by best case scenario)
      if (min_score==INT_MAX && current_score+pattern_left[top_level] > text_left+max_distance) return INT_MAX;
    }
  }
  // Return results
  return min_score;
}






