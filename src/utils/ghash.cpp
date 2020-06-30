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
 * DESCRIPTION: Simple Hash Implementation (Generic Key, Generic Value)
 */

#include "../utils/ghash.h"

/*
 * Constants
 */
#define GHASH_SIZE_PER_ELEMENT 56

/*
 * Allocators
 */
#undef  uthash_malloc
#define uthash_malloc(sz)   (ghash->mm_allocator!=NULL ? mm_allocator_malloc(ghash->mm_allocator,sz) : malloc(sz))
#undef  uthash_free
#define uthash_free(ptr,sz) (ghash->mm_allocator!=NULL) ? mm_allocator_free(ghash->mm_allocator,ptr) : free(ptr);

/*
 * Constructor
 */
ghash_t* ghash_new(
    const int key_size,
    mm_allocator_t* const mm_allocator) {
  ghash_t* const ghash =  (mm_allocator==NULL) ?
      (ghash_t*) malloc(sizeof(ghash_t)) :
       mm_allocator_alloc(mm_allocator, ghash_t);
  ghash->head = NULL; // uthash initializer
  ghash->key_size = key_size;
  ghash->mm_allocator = mm_allocator;
  return ghash;
}
void ghash_clear(ghash_t* const ghash) {
  ghash_element_t *ghash_element, *tmp;
  HASH_ITER(hh,ghash->head,ghash_element,tmp) {
    HASH_DEL(ghash->head,ghash_element);
    uthash_free(ghash_element->key,0);
    uthash_free(ghash_element,0);
  }
  ghash->head = NULL; // uthash initializer
}
void ghash_delete(ghash_t* const ghash) {
  ghash_clear(ghash);
  uthash_free(ghash,0);
}
/*
 * Handlers (Type-unsafe Accessors)
 */
ghash_element_t* ghash_handler_get_element(
    ghash_t* const ghash,
    const void* const key) {
  ghash_element_t* ghash_element;
  HASH_FIND(hh,ghash->head,key,ghash->key_size,ghash_element);
  return ghash_element;
}
void* ghash_handler_get(
    ghash_t* const ghash,
    const void* const key) {
  ghash_element_t* const ghash_element = ghash_handler_get_element(ghash,key);
  return (ghash_element!=NULL) ? ghash_element->element : NULL;
}
void ghash_handler_insert(
    ghash_t* ghash,
    const void* const key,
    void* const element) {
  ghash_element_t* ghash_element = ghash_handler_get_element(ghash,key);
  if (ghash_element==NULL) {
    // Allocate
    ghash_element = (ghash_element_t*) uthash_malloc(sizeof(ghash_element_t));
    // Set key/element
    const int key_size = ghash->key_size;
    ghash_element->key = uthash_malloc(key_size);
    memcpy(ghash_element->key,key,key_size);
    ghash_element->element = element;
    // Add to hash
    HASH_ADD_KEYPTR(hh,ghash->head,ghash_element->key,key_size,ghash_element);
  } else {
    // Replace element (no free whatsoever)
    ghash_element->element = element;
  }
}
void ghash_handler_remove(
    ghash_t* ghash,
    const void* const key) {
  ghash_element_t* ghash_element = ghash_handler_get_element(ghash,key);
  if (ghash_element) {
    HASH_DEL(ghash->head,ghash_element);
    uthash_free(ghash_element->key,0);
    uthash_free(ghash_element,0);
  }
}
/*
 * Accessors
 */
bool ghash_is_contained(
    ghash_t* const ghash,
    const void* const key) {
  return (ghash_handler_get_element(ghash,key)!=NULL);
}
int ghash_get_num_elements(ghash_t* const ghash) {
  return (int)HASH_COUNT(ghash->head);
}
int ghash_get_size(ghash_t* const ghash) {
  /*
   * The hash handle consumes about 32 bytes per item on a 32-bit system, or 56 bytes per item on a 64-bit system.
   * The other overhead costs (the buckets and the table) are negligible in comparison.
   * You can use HASH_OVERHEAD to get the overhead size, in bytes, for a hash table.
   */
  return ghash_get_num_elements(ghash)*GHASH_SIZE_PER_ELEMENT;
}
/*
 * Iterator
 */
ghash_iterator_t* ghash_iterator_new(ghash_t* const ghash) {
  // Allocate
  ghash_iterator_t* const iterator = (ghash_iterator_t*) uthash_malloc(sizeof(ghash_iterator_t));
  // Init
  iterator->ghash = ghash;
  iterator->current = NULL;
  iterator->next = ghash->head;
  return iterator;
}
void ghash_iterator_delete(ghash_iterator_t* const iterator) {
  ghash_t* const ghash = iterator->ghash;
  uthash_free(iterator,0);
}
bool ghash_iterator_eoi(ghash_iterator_t* const iterator) {
  return (iterator->next != NULL);
}
bool ghash_iterator_next(ghash_iterator_t* const iterator) {
  if (iterator->next != NULL) {
    iterator->current = iterator->next;
    iterator->next = (ghash_element_t*) iterator->next->hh.next;
    return true;
  } else {
    iterator->current = NULL;
    return false;
  }
}
void* ghash_iterator_get_key(ghash_iterator_t* const iterator) {
  return iterator->current->key;
}
void* ghash_iterator_get_element(ghash_iterator_t* const iterator) {
  return iterator->current->element;
}

