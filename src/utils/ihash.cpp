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
 * DESCRIPTION: Simple Hash Implementation (Interger Key, Generic Value)
 */

#include "../utils/ihash.h"

/*
 * Constants
 */
#define IHASH_SIZE_PER_ELEMENT 56

/*
 * Allocators
 */
#undef  uthash_malloc
#define uthash_malloc(sz) (ihash->mm_allocator!=NULL ? mm_allocator_malloc(ihash->mm_allocator,sz) : malloc(sz))
#undef  uthash_free
#define uthash_free(ptr,sz) (ihash->mm_allocator!=NULL) ? mm_allocator_free(ihash->mm_allocator,ptr) : free(ptr);

/*
 * Constructor
 */
ihash_t* ihash_new(mm_allocator_t* const mm_allocator) {
  ihash_t* const ihash = (mm_allocator==NULL) ?
      (ihash_t*) malloc(sizeof(ihash_t)) :
      mm_allocator_alloc(mm_allocator,ihash_t);
  ihash->head = NULL; // uthash initializer
  ihash->mm_allocator = mm_allocator;
  return ihash;
}
void ihash_clear(ihash_t* const ihash) {
  ihash_element_t *ihash_element, *tmp;
  HASH_ITER(hh,ihash->head,ihash_element,tmp) {
    HASH_DEL(ihash->head,ihash_element);
    uthash_free(ihash_element,0);
  }
  ihash->head = NULL; // uthash initializer
}
void ihash_delete(ihash_t* const ihash) {
  ihash_clear(ihash);
  uthash_free(ihash,0);
}
/*
 * Handlers (Type-unsafe Accessors)
 */
ihash_element_t* ihash_handler_get_element(
    ihash_t* const ihash,
    const int key) {
  ihash_element_t *ihash_element;
  HASH_FIND_INT(ihash->head,&key,ihash_element);
  return ihash_element;
}
void* ihash_handler_get(
    ihash_t* const ihash,
    const int key) {
  ihash_element_t* const ihash_element = ihash_handler_get_element(ihash,key);
  return (ihash_element!=NULL) ? ihash_element->element : NULL;
}
void ihash_handler_insert(
    ihash_t* ihash,
    const int key,
    void* const element) {
  ihash_element_t* ihash_element = ihash_handler_get_element(ihash,key);
  if (ihash_element==NULL) {
    // Allocate
    ihash_element = (ihash_element_t*) uthash_malloc(sizeof(ihash_element_t));
    // Set key/element
    ihash_element->key = key;
    ihash_element->element = element;
    // Add to hash
    HASH_ADD_INT(ihash->head,key,ihash_element);
  } else {
    // Replace element (no free whatsoever)
    ihash_element->element = element;
  }
}
void ihash_handler_remove(
    ihash_t* ihash,
    const int key) {
  ihash_element_t* ihash_element = ihash_handler_get_element(ihash,key);
  if (ihash_element) {
    HASH_DEL(ihash->head,ihash_element);
    uthash_free(ihash_element,0);
  }
}
/*
 * Accessors
 */
bool ihash_is_contained(
    ihash_t* const ihash,
    const int key) {
  return (ihash_handler_get_element(ihash,key)!=NULL);
}
int ihash_get_num_elements(ihash_t* const ihash) {
  return (int)HASH_COUNT(ihash->head);
}
int ihash_get_size(ihash_t* const ihash) {
  /*
   * The hash handle consumes about 32 bytes per item on a 32-bit system, or 56 bytes per item on a 64-bit system.
   * The other overhead costs (the buckets and the table) are negligible in comparison.
   * You can use HASH_OVERHEAD to get the overhead size, in bytes, for a hash table.
   */
  return ihash_get_num_elements(ihash)*IHASH_SIZE_PER_ELEMENT;
}
/*
 * Miscellaneous
 */
int ihash_cmp_keys(int* a,int* b) {
  /*
   * return (int) -1 if (a < b)
   * return (int)  0 if (a == b)
   * return (int)  1 if (a > b)
   */
  return *a-*b;
}
#define ihash_cmp_keys_wrapper(arg1,arg2) ihash_cmp_keys((int*)arg1,(int*)arg2)
void ihash_sort_by_key(ihash_t* ihash) {
  HASH_SORT(ihash->head,ihash_cmp_keys_wrapper); // Sort
}
/*
 * Iterator
 */
ihash_iterator_t* ihash_iterator_new(ihash_t* const ihash) {
  // Allocate
  ihash_iterator_t* const iterator = (ihash_iterator_t*) uthash_malloc(sizeof(ihash_iterator_t));
  // Init
  iterator->ihash = ihash;
  iterator->current = NULL;
  iterator->next = ihash->head;
  return iterator;
}
void ihash_iterator_delete(ihash_iterator_t* const iterator) {
  ihash_t* const ihash = iterator->ihash;
  uthash_free(iterator,0);
}
bool ihash_iterator_eoi(ihash_iterator_t* const iterator) {
  return (iterator->next != NULL);
}
bool ihash_iterator_next(ihash_iterator_t* const iterator) {
  if (iterator->next != NULL) {
    iterator->current = iterator->next;
    iterator->next = (ihash_element_t*) iterator->next->hh.next;
    return true;
  } else {
    iterator->current = NULL;
    return false;
  }
}
int ihash_iterator_get_key(ihash_iterator_t* const iterator) {
  return iterator->current->key;
}
void* ihash_iterator_get_element(ihash_iterator_t* const iterator) {
  return iterator->current->element;
}
