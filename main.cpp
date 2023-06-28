#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <tuple>
#include <vector>

std::string toBinaryString(uint64_t n) {
  int num_bits = sizeof(uint64_t) * 8; //initializes a variable num_bits to hold the number of bits in a uint64_t, bytes to bits
  std::string str;
  for (int i = 0; i < num_bits; i++) {
    str += ((n & 1) + '0');   //bitwise 'and' operation with 1 , then rightshift
    n >>= 1;          // n is the least significant bit
  }
  std::reverse(str.begin(), str.end()); //reverse string
  return str;
}
//create a bitmask with a range of set bits
uint64_t create_mask(uint64_t a, uint64_t b) { //, a and b, which represent the starting and ending positions of the set bits in the mask.
  uint64_t r = 0;
  for (unsigned i = a; i <= b; i++) {
    r |= 1 << i; //bitwise or
  }
  return r;
}

uint64_t extract_bits(uint64_t num, uint64_t start, uint64_t stop) {
  return (num & create_mask(start, stop)) >> start; // clears all the bits outside desored range, while preserving the bits within the range.
}                                                   //then right shifted by start to align extracted bits to rightmost position

bool IsPowerOfTwo(int x) { return (x != 0) && ((x & (x - 1)) == 0); }

int main(int argc, char **argv) {
  // argc represents the number of command-line arguments
  // argv is an array of pointers to command-line argument strings
  // checks that all parameters are given
  if (argc < 7) {
    std::cout << argv[0]
              << " BLOCKSIZE L1_SIZE L1_ASSOC L2_SIZE L2_ASSOC "
                 "memory_trace_files/trace1.txt"
              << std::endl;
    return EXIT_FAILURE;
  }
  // extract the parameters and convert them from str to int
  const int BLOCKSIZE = std::atoi(argv[1]);
  const int L1_SIZE = std::atoi(argv[2]);
  const int L1_ASSOC = std::atoi(argv[3]);
  const int L2_SIZE = std::atoi(argv[4]);
  const int L2_ASSOC = std::atoi(argv[5]);

  // make sure all the parameters are power of two
  if (!IsPowerOfTwo(BLOCKSIZE)) {
    std::cerr << "BLOCKSIZE: " << BLOCKSIZE << " must be power of two"
              << std::endl;
    return EXIT_FAILURE;
  }
  if (!IsPowerOfTwo(L1_SIZE)) {
    std::cerr << "L1_SIZE: " << L1_SIZE << " must be power of two" << std::endl;
    return EXIT_FAILURE;
  }
  if (!IsPowerOfTwo(L1_ASSOC)) {
    std::cerr << "L1_ASSOC: " << L1_ASSOC << " must be power of two"
              << std::endl;
    return EXIT_FAILURE;
  }
  if (!IsPowerOfTwo(L2_SIZE)) {
    std::cerr << "L2_SIZE: " << L2_SIZE << " must be power of two" << std::endl;
    return EXIT_FAILURE;
  }
  if (!IsPowerOfTwo(L2_ASSOC)) {
    std::cerr << "L2_ASSOC: " << L2_ASSOC << " must be power of two"
              << std::endl;
    return EXIT_FAILURE;
  }

  // open given file
  std::ifstream infile{argv[6]}; //infile is an instance used to read data
  if (!infile) {
    std::cerr << "file: " << argv[6] << " could not be opened" << std::endl;
    return EXIT_FAILURE;
  }

  const int OFFSET_BITS = std::log2(BLOCKSIZE);
  const int L1_NUM_SETS = L1_SIZE / (BLOCKSIZE * L1_ASSOC);
  const int L1_INDEX_BITS = std::log2(L1_NUM_SETS);
  
  //make 2-D L1 cache
  std::vector<std::vector<std::tuple<uint64_t, uint64_t, bool, uint64_t>>>
      L1_cache(L1_NUM_SETS);              // tag, recent use score, dirty bit, and address
  for (int i = 0; i < L1_NUM_SETS; ++i) { //initialise all elements to (0, 0, false, 0)
    for (int j = 0; j < L1_ASSOC; ++j) {
      L1_cache[i].push_back(std::make_tuple(0, 0, false, 0));
    }
  }
  // L1 counters
  int L1_READS = 0;
  int L1_READ_MISSES = 0;
  int L1_WRITES = 0;
  int L1_WRITE_MISSES = 0;
  double L1_MISS_RATE = 0;
  int L1_WRITE_BACKS = 0;

  const int L2_NUM_SETS = L2_SIZE / (BLOCKSIZE * L2_ASSOC);
  const int L2_INDEX_BITS = std::log2(L2_NUM_SETS);

  //make 2-D L2 cache
  std::vector<std::vector<std::tuple<uint64_t, uint64_t, bool, uint64_t>>>
      L2_cache(L2_NUM_SETS);
  for (int i = 0; i < L2_NUM_SETS; ++i) {
    for (int j = 0; j < L2_ASSOC; ++j) {
      L2_cache[i].push_back(std::make_tuple(0, 0, false, 0));
    }
  }
  //L2 counters
  int L2_READS = 0;
  int L2_READ_MISSES = 0;
  int L2_WRITES = 0;
  int L2_WRITE_MISSES = 0;
  double L2_MISS_RATE = 0;
  int L2_WRITE_BACKS = 0;

  std::string line;
  int line_number = 0;
  // read file line by line
  while (std::getline(infile, line)) {
    // to keep track of recently used.
    line_number++;
    std::istringstream iss(line);//initialise line as a string to extract values
    char instruction;
    uint64_t address;

    // read get params from the line
    if (!(iss >> instruction >> std::hex >> address)) { //stream interpretation, address variable will be assigned the integer value corresponding to the hexadecimal input 
      std::cerr << "there was an error while parsing file" << std::endl;
      break;
    }
    //get index and tag bits
    uint64_t l1_index_bits =
        extract_bits(address, OFFSET_BITS, OFFSET_BITS + L1_INDEX_BITS - 1);
    uint64_t l1_tag_bits = extract_bits(address, L1_INDEX_BITS, 63);

    switch (instruction) {
    case 'r': {
      // std::cout << "r" << std::endl;
      // search L1 cache
      auto found_in_L1 = std::find_if(
          L1_cache[l1_index_bits].begin(), L1_cache[l1_index_bits].end(),
          [&l1_tag_bits](const auto &tuple) {
            return std::get<0>(tuple) == l1_tag_bits;
          });
      if (found_in_L1 != L1_cache[l1_index_bits].end()) { //condition verifies if a matching element was found in the L1 cache.
        // found in L1
        L1_READS++;
        // update its recent use score
        std::get<1>(*found_in_L1) = line_number;
      } else {
        // NOT found in L1
        L1_READ_MISSES++;
        // something must now be evicted
        // find least recently used from L1 to evict
        // with lowest line_number
        auto least_recenlty_used_in_l1_iter = std::min_element(
            L1_cache[l1_index_bits].begin(), L1_cache[l1_index_bits].end(),
            [](const auto &tuple_a, const auto &tuple_b) {
              return std::get<1>(tuple_a) < std::get<1>(tuple_b);
            });

        // save a copy for L2
        auto evicted_l1 = *least_recenlty_used_in_l1_iter;

        // search L2 cache
        uint64_t l2_index_bits =
            extract_bits(address, OFFSET_BITS, OFFSET_BITS + L2_INDEX_BITS - 1);
        uint64_t l2_tag_bits = extract_bits(address, L2_INDEX_BITS, 63);

        auto found_in_L2 = std::find_if(
            L2_cache[l2_index_bits].begin(), L2_cache[l2_index_bits].end(),
            [&l2_tag_bits](const auto &tuple) {
              return std::get<0>(tuple) == l2_tag_bits;
            });
        if (found_in_L2 != L2_cache[l2_index_bits].end()) {
          // found in L2
          L2_READS++;
          // was it dirty ? update dirty bit of l1 block
          std::get<2>(*least_recenlty_used_in_l1_iter) =
              std::get<2>(*found_in_L2);

        } else {
          // NOT found in L2 as well
          L2_READ_MISSES++;
          // if L1 evicted is dirty write back
          if (std::get<2>(*least_recenlty_used_in_l1_iter)) {
            L1_WRITE_BACKS++;
          }
          // definetely not dirty
          std::get<2>(*least_recenlty_used_in_l1_iter) = false;
        }

        // finish eviction from l1
        // mind that tag bits for the caches are not the same
        std::get<0>(*least_recenlty_used_in_l1_iter) = l1_tag_bits;
        std::get<1>(*least_recenlty_used_in_l1_iter) = line_number;
        std::get<3>(*least_recenlty_used_in_l1_iter) = address;

        // now evict least recently used from L2
        auto least_recenlty_used_l2 = std::min_element(
            L2_cache[l2_index_bits].begin(), L2_cache[l2_index_bits].end(),
            [](const auto &tuple_a, const auto &tuple_b) {
              return std::get<1>(tuple_a) < std::get<1>(tuple_b);
            });

        // if L2 evicted is dirty write back
        if (std::get<2>(*least_recenlty_used_l2)) {
          L2_WRITE_BACKS++;
        }
        // now replace it with memory evicted from L1
        // we have to recalculate address of just evicted block from L1
        std::get<0>(*least_recenlty_used_l2) =
            extract_bits(std::get<3>(evicted_l1), L2_INDEX_BITS, 63);
        std::get<1>(*least_recenlty_used_l2) = line_number;
        std::get<2>(*least_recenlty_used_l2) = std::get<2>(evicted_l1);
      }
    } break;
    case 'w': {
      // std::cout << "r" << std::endl;
      // search L1 cache
      auto found_in_L1 = std::find_if(
          L1_cache[l1_index_bits].begin(), L1_cache[l1_index_bits].end(),
          [&l1_tag_bits](const auto &tuple) {
            return std::get<0>(tuple) == l1_tag_bits;
          });
      if (found_in_L1 != L1_cache[l1_index_bits].end()) {
        // found in L1
        L1_WRITES++;
        // update its recent use score
        std::get<1>(*found_in_L1) = line_number;
        // it is now dirty
        std::get<2>(*found_in_L1) = true;
      } else {
        // NOT found in L1
        L1_WRITE_MISSES++;
        // evict and find least recently used with lowest line_number
        auto least_recenlty_used_l1_iter = std::min_element(
            L1_cache[l1_index_bits].begin(), L1_cache[l1_index_bits].end(),
            [](const auto &tuple_a, const auto &tuple_b) {
              return std::get<1>(tuple_a) < std::get<1>(tuple_b);
            });
        // save a copy for L2
        auto evicted_l1 = *least_recenlty_used_l1_iter;

        // search L2 cache
        uint64_t l2_index_bits =
            extract_bits(address, OFFSET_BITS, OFFSET_BITS + L2_INDEX_BITS - 1);
        uint64_t l2_tag_bits = extract_bits(address, L2_INDEX_BITS, 63);

        auto found_in_L2 = std::find_if(
            L2_cache[l2_index_bits].begin(), L2_cache[l2_index_bits].end(),
            [&l2_tag_bits](const auto &tuple) {
              return std::get<0>(tuple) == l2_tag_bits;
            });
        if (found_in_L2 != L2_cache[l2_index_bits].end()) {
          // found in L2
          L2_WRITES++;
          // now add from L2  to the L1 cache
          // mind that tag bits for the caches are not the same

          // if L2 evicted is dirty write back
          if (std::get<2>(*found_in_L2)) {
            L2_WRITE_BACKS++;
          }
          // was l1 dirty?
          std::get<2>(*found_in_L2) = std::get<2>(evicted_l1);
        } else {
          // NOT found in L2 either
          L2_WRITE_MISSES++;
          // if about to be evicted from L1 was dirty
          if (std::get<2>(*least_recenlty_used_l1_iter)) {
            L1_WRITE_BACKS++;
          }

          std::get<2>(*least_recenlty_used_l1_iter) = false;
        }
        // finish eviction from l1
        // mind that tag bits for the caches are not the same
        std::get<0>(*least_recenlty_used_l1_iter) = l1_tag_bits;
        std::get<1>(*least_recenlty_used_l1_iter) = line_number;
        std::get<3>(*least_recenlty_used_l1_iter) = address;

        // now evict least recently used from L2
        auto least_recenlty_used_l2 = std::min_element(
            L2_cache[l2_index_bits].begin(), L2_cache[l2_index_bits].end(),
            [](const auto &tuple_a, const auto &tuple_b) {
              return std::get<1>(tuple_a) < std::get<1>(tuple_b);
            });

        // if L2 evicted is dirty write back
        if (std::get<2>(*least_recenlty_used_l2)) {
          L2_WRITE_BACKS++;
        }
        // now replace it with memory evicted from L1
        // we have to recalculate address of just evicted from L1
        std::get<0>(*least_recenlty_used_l2) =
            extract_bits(std::get<3>(evicted_l1), L2_INDEX_BITS, 63);
        std::get<1>(*least_recenlty_used_l2) = line_number;
        std::get<2>(*least_recenlty_used_l2) = std::get<2>(evicted_l1);
      }
      break;
    }
    default:
      break;
    }
  }
  // report
  L1_MISS_RATE = (double)(L1_READ_MISSES + L1_WRITE_MISSES) / line_number;
  L2_MISS_RATE = (double)(L2_READ_MISSES + L2_WRITE_MISSES) /
                 (L1_READ_MISSES + L1_WRITE_MISSES);
  std::cout << "i.    number of L1 reads: " << L1_READS << std::endl;
  std::cout << "ii.   number of L1 read misses: " << L1_READ_MISSES
            << std::endl;
  std::cout << "iii.  number of L1 writes: " << L1_WRITES << std::endl;
  std::cout << "iv.   number of L1 write misses: " << L1_WRITE_MISSES
            << std::endl;
  std::cout << "v.    L1 miss rate: " << L1_MISS_RATE << std::endl;
  std::cout << "vi.   number of writebacks from L1 memory: " << L1_WRITE_BACKS
            << std::endl;

  std::cout << "vii.  number of L2 reads: " << L2_READS << std::endl;
  std::cout << "viii. number of L2 read misses: " << L2_READ_MISSES
            << std::endl;
  std::cout << "ix.   number of L2 writes:" << L2_WRITES << std::endl;
  std::cout << "x.    number of L2 write misses: " << L2_WRITE_MISSES
            << std::endl;
  std::cout << "xi.   L2 miss rate: " << L2_MISS_RATE << std::endl;
  std::cout << "xii.  number of writebacks from L2 memory: " << L2_WRITE_BACKS
            << std::endl;

  // int total_time = ((L1_READS + L1_WRITES) * 1) +
  //                  ((L2_READS + L2_WRITES) * 20) +
  //                  ((L2_READ_MISSES + L2_WRITE_MISSES) * 200);
  // std::cout << "total access time: " << total_time << " ns" << std::endl;

  return EXIT_SUCCESS;
}
