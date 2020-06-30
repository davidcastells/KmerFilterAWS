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
 * DESCRIPTION:
 *   Filter based on general k-mer counting as to quickly filter out
 *   candidates that cannot align against its region-text
 */

#include "kmer_filter.h"

#include "../filter/pattern_tiling.h"
#include "../utils/dna_text.h"

/*
 * Constants                                        Histograms size
 */
#define KMER_COUNTING_MASK_3   0x000000000000003Full  /* 256 B */
#define KMER_COUNTING_MASK_4   0x00000000000000FFull  /*  1 KB */
#define KMER_COUNTING_MASK_5   0x00000000000003FFull  /*  4 KB */
#define KMER_COUNTING_MASK_6   0x0000000000000FFFull  /* 16 KB */
#define KMER_COUNTING_MASK_7   0x0000000000003FFFull  /* 64 KB */
#define KMER_COUNTING_MASK_8   0x000000000000FFFFull  /* 256 KB */
#define KMER_COUNTING_MASK_9   0x000000000003FFFFull  /*   1 MB */
#define KMER_COUNTING_MASK_10  0x00000000000FFFFFull  /*   4 MB */
#define KMER_COUNTING_MASK_11  0x00000000003FFFFFull  /*  16 MB */
#define KMER_COUNTING_MASK_12  0x0000000000FFFFFFull  /*  64 MB */
#define KMER_COUNTING_MASK_13  0x0000000003FFFFFFull  /* 256 MB */

#define KMER_COUNTING_MASK_INDEX(kmer_idx) \
  ((kmer_idx) & kmer_counting->kmer_mask)
#define KMER_COUNTING_ADD_INDEX(kmer_idx,character) \
  kmer_idx = (kmer_idx<<2 | (enc_char))
#define KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char) \
  kmer_idx = KMER_COUNTING_MASK_INDEX(kmer_idx<<2 | (enc_char))

/*
 * Uncalled bases handling
 */
#define KMER_COUNTING_FILTER_CHAR(character) (dna_encode(character) % ENC_DNA_CHAR_N) // FIXME: Correct-lossless?

/*
 * Setup
 */
kmer_counting_nway_t* kmer_counting_new(const uint64_t kmer_length, mm_allocator_t* const mm_allocator) 
{
  // Allocate
  kmer_counting_nway_t* const kmer_counting = mm_allocator_alloc(mm_allocator,kmer_counting_nway_t);
  // Filter parameters
  kmer_counting->kmer_length = kmer_length;
  switch (kmer_length) { // Check kmer length
    case 3: kmer_counting->kmer_mask = KMER_COUNTING_MASK_3; break;
    case 4: kmer_counting->kmer_mask = KMER_COUNTING_MASK_4; break;
    case 5: kmer_counting->kmer_mask = KMER_COUNTING_MASK_5; break;
    case 6: kmer_counting->kmer_mask = KMER_COUNTING_MASK_6; break;
    case 7: kmer_counting->kmer_mask = KMER_COUNTING_MASK_7; break;
    case 8: kmer_counting->kmer_mask = KMER_COUNTING_MASK_8; break;
    case 9: kmer_counting->kmer_mask = KMER_COUNTING_MASK_9; break;
    case 10: kmer_counting->kmer_mask = KMER_COUNTING_MASK_10; break;
    case 11: kmer_counting->kmer_mask = KMER_COUNTING_MASK_11; break;
    case 12: kmer_counting->kmer_mask = KMER_COUNTING_MASK_12; break;
    case 13: kmer_counting->kmer_mask = KMER_COUNTING_MASK_13; break;
    default:
      fprintf(stderr,"K-mer counting. Invalid proposed k-mer length\n");
      exit(1);
      break;
  }
  kmer_counting->num_kmers = POW4(kmer_counting->kmer_length);
  // Allocate histogram tables
  const uint64_t kmer_table_size = kmer_counting->num_kmers * sizeof(kmer_count_int_t);
  void* const memory = mm_allocator_calloc(mm_allocator,2*kmer_table_size,uint8_t,true);
  kmer_counting->kmer_count_text = (kmer_count_int_t*) memory;
  kmer_counting->kmer_count_pattern = (kmer_count_int_t*)(memory + kmer_table_size);
  // MM
  kmer_counting->mm_allocator = mm_allocator;
  // Return
  return kmer_counting;
}
void kmer_counting_destroy(
    kmer_counting_nway_t* const kmer_counting) {
  mm_allocator_free(kmer_counting->mm_allocator,kmer_counting->kmer_count_text);
  mm_allocator_free(kmer_counting->mm_allocator,kmer_counting);
}
/*
 * Pattern prepare
 * @param kmer_counting
 * @param key
 * @param key_length    is the number of bases of the input sequence
 * 
 * We do an sliding window over the key. The n-gram is stored in kmer_idx.
 * For every incoming character we shift left and add the new symbol
 */
void kmer_counting_pattern_compute_histogram(kmer_counting_nway_t* const kmer_counting,
    uint8_t* const key,
    const uint64_t key_length) 
{
  // Parameters
  kmer_count_int_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  
  // Set key parameters
  kmer_counting->key = key;
  kmer_counting->key_length = key_length;
  kmer_counting->num_key_kmers = key_length - (kmer_counting->kmer_length-1);
  
  // Count until chunk end
  uint64_t pos, kmer_idx = 0, acc = 0;
  
  for (pos=0;pos<key_length;++pos) 
  {
    const uint8_t character = key[pos];
    if (character == DNA_CHAR_N) 
    {
      acc = 0;
    }
    else 
    {
      KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,dna_encode(character)); // Update kmer-index
      
      if (acc < kmer_counting->kmer_length-1) 
      {
        ++acc; // Inc accumulator
      } 
      else 
      {
        // Increment kmer-count (kmer-counters for all tiles are stored together)
        ++(kmer_count_pattern[kmer_idx]);
      }
    }
  }
}


/**
 * K-mer counting
 */
uint64_t kmer_counting_min_bound(
    kmer_counting_nway_t* const kmer_counting,
    const uint8_t* const text,
    const uint64_t text_length,
    const uint64_t max_error) {
  // Parameters
  const uint64_t kmer_length = kmer_counting->kmer_length;
  kmer_count_int_t* const kmer_count_pattern = kmer_counting->kmer_count_pattern;
  kmer_count_int_t* const kmer_count_text = kmer_counting->kmer_count_text;
  uint64_t kmer_idx = 0, kmer_end, kmer_begin;
  // Prepare filter
  memset(kmer_counting->kmer_count_text,0,kmer_counting->num_kmers*sizeof(kmer_count_int_t));
  // Prepare text
  kmer_counting->max_text_kmers = 0;
  kmer_counting->curr_text_kmers = 0;
  // Initial fill (kmer)
  for (kmer_end=0;kmer_end<kmer_length-1;++kmer_end) {
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,KMER_COUNTING_FILTER_CHAR(text[kmer_end]));
  }
  // Sliding window
  for (kmer_begin=0;kmer_end<text_length;++kmer_begin,++kmer_end) {
    // Fetch counters & store them in window
    const uint8_t enc_char = KMER_COUNTING_FILTER_CHAR(text[kmer_end]);
    KMER_COUNTING_ADD_INDEX__MASK(kmer_idx,enc_char);
    const uint64_t kmer_offset = kmer_idx;
    kmer_count_int_t* const text_count_ptr = kmer_count_text + kmer_offset;
    kmer_count_int_t* const pattern_count_ptr = kmer_count_pattern + kmer_offset;
    // Increment kmer counts
    const kmer_count_int_t text_count = *text_count_ptr;
    const kmer_count_int_t pattern_count = *pattern_count_ptr;
    if (pattern_count > 0 && text_count < pattern_count) {
      ++(kmer_counting->curr_text_kmers);
      kmer_counting->max_text_kmers = MAX(kmer_counting->max_text_kmers,kmer_counting->curr_text_kmers);
    }
    ++(*text_count_ptr);
  }
  // Compute min-error bound
  const uint64_t kmer_diff = kmer_counting->num_key_kmers - kmer_counting->max_text_kmers;
  return DIV_CEIL(kmer_diff,kmer_length);
}


