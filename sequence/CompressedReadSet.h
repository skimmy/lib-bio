#ifndef COMPRESSED_READ_SET_H
#define COMPRESSED_READ_SET_H

#include "CompressedSequence.h"
#include "Read.h"

#include <cstddef>
#include <vector>

// TODO LIST
// 1. Perform resizing of the read size in order to have a multiple of 8 sequence size (in bits) -- DONE
// 2. Add the 'append' methods (and the += operator as well)
//    append takes Read and vector<Read> as parameter (i.e. two overload of the method) -- DONE
// 3. Add 'loadFromFile' methods --  DONE
// 4. Add 'getReads(begin, count)' method to retrieve the 'count' reads [begin, begin+count-1] -- DONE


/**
 * \brief This class represents a compressed set of reads. 
 *
 * A set of reads is a sequence of homogeneous reads (i.e. same 
 * length). The compression is performed using the DNACompressedSymbol
 * and the number of bits for each compressed symbol is given by
 * the static variable CompressedReadSet::BaseBitLength (default value is 4).
 * For efficeincy reasons read are eventually padded to reach
 * the first value x such that \f$ BaseBitLengh*x \f$ is a multiple if
 * 8, this is due to internal representation of the compressed
 * base which are 8 bits words, the pad is filled with the '0'
 * value which is the sequence made of BaseBitLength zeros.
 */
class CompressedReadSet : public CompressedSequence {
 public:
  /**
   * \brief The number of bits per base (default value is 4)
   */
  static int BaseBitLength;
  /**
   * \brief The length of the reads when used by the default constructor
   */
  static int DefaultReadSize;

 private:
  size_t readCount;
  size_t readLength;
 public:
  // ---------------------------------------------------------
  // -------------- CONSTRUCTORS AND DESTRUCTOR --------------
  // ---------------------------------------------------------

  /**
   * \brief Creates an emtpy set. 
   *
   * The interenal sequence is therefore empty (yet initialized). The read
   * size is initialized using DefaultReadSize static member.
   * Note that read size could not be changed afterwards, you
   * should consider using the constructor allowing to specify
   * the size of the read.
   */
  CompressedReadSet();
  /**
   * \brief Creates an empty set where the size of a single
   * read is specified
   *
   * This constructor works exactly as the default one (i.e.
   * creates an emtpy sequence of 0 reads) but it requires the
   * programmer to specify the size of a read (rather than using
   * the default value).
   *
   * \param readSize The size of the reads
   *
   * \sa CompressedReadSet()
   */
  CompressedReadSet(size_t readSize);
  /**
   * \brief  Constructs the CompressedReadSet from an array of Read.
   *
   * The number of reads is determined by the parameter given
   * and the length of a single read is obtained from the read
   * of index 0 of the passed array. It is therefore important
   * that all reads in the input array are homogenous (all have
   * the same length) as prescribed by the class.
   *
   * \param reads the array of Read objects
   * \param n the number of Read objects in reads
   *
   * \sa CompressedReadSet(const vector<Read>& reads)
   */
  CompressedReadSet(const Read* reads, size_t n);
  /**
   * \brief Constructs the CompressedReadSet form a vector<Read> object.
   * 
   * The number of reads is obtained by the size of the input
   * vector while the length of the single read is obtained
   * by the 0 index read of the input vector. It is therefore important
   * that all reads in the input array are homogenous (all have
   * the same length) as prescribed by the class.
   *
   * \param reads the vector of reads
   *
   * \sa CompressedReadSet(const Read* reads, size_t n)
   */
  CompressedReadSet(const vector<Read>& reads);
  /**
   * \brief Copy constructor, creates a copy of the input CompressedReadSet
   *
   * \param other The CompressedReadSet to be copied
   */
  CompressedReadSet(const CompressedReadSet& other);
  ~CompressedReadSet();

  // ---------------------------------------------------------
  // ------------------ GET AND SET METHODS ------------------
  // ---------------------------------------------------------

  /**
   * \brief Returns the number of read stored in the set.
   *
   * \return The number of reads stored in the set
   */
  size_t getReadCount() const;

  /**
   * \brief Returns a copy of the sequence consisting of count reads
   * from the position start. 
   *
   * If start + count - 1 is greater
   * then the actual number of reads the remaining reads are
   * filled with all zeros. The returned sequence is synamically
   * allocate array, programmers should take care of freeing
   * associated memory when sequence is no more needed.
   *   
   * \param start The index of first read to be returned
   * \param count The total number of reads returned
   * \return a copy of the sequence for returned reads
   * possibily padded with all zeros reads
   *
   */
  uint8_t* getReads(size_t start, size_t count) const;

  // ---------------------------------------------------------
  //                      MODIFY METHODS 
  // ---------------------------------------------------------
  /**
   * \brief Appends a new read to the sequence.
   *
   * The new read should be heomogeneous in length with all the
   * already inserted reads. This means that when appending the
   * read its size is translated to match the actual size. If
   * the real size of read is smaller zeros are padded at the end
   * otherwise base sequence is truncated to actual length.
   * Since the compressed seqeunce is enlarged for the amount of
   * space required by the new read, it is more convenient to use
   * other overloaded versions of this method when more reads need
   * to be appended to the compressed sequence.
   *
   * \param read The read to be appended
   * \return This object after append has been performed
   *
   * \sa CompressedReadSet& append(const vector<Read>& reads)
   *
   */
  CompressedReadSet& append(const Read& read);

  /**
   * \brief Appends several reads to the sequence.
   *
   * The new reads should be heomogeneous in length with all the
   * already inserted reads. This means that when appending them
   * their size is translated to match the actual size. If
   * the real size of such reads is smaller zeros are padded at the end
   * otherwise base sequences are truncated to actual length.
   *
   * \param reads The vector containing the reads to be appended
   * \return This object after append has been performed
   *
   * \sa CompressedReadSet& append(const Read& read)
   */
  CompressedReadSet& append(const vector<Read>& reads);

  // I/O methods
  /**
   * \brief Writes the content to a file. 
   *
   * In the target file are stored
   * *in binary format) the size of the entire sequence, the number
   * of reads and the length of each read, after this information
   * the raw sequence is stored (again in binary format).
   * 
   * \param fileName the full path of the file.
   *
   * \sa CompressedReadSet& loadFromFile(const string& fileName)
   */
  void writeToFile(const string& fileName) const;

  /**
   * \brief Loads the content of fileNmae. 
   *
   * The target file must follow
   * the rules defined in writeToFile methods. The content of
   * this object is destroyed after loading new content from
   * file and is no longer available when the method returns.
   *
   * \param fileName The full path of the file
   * \return This object after loading has been performed
   *
   * \sa writeToFile(const string& fileName) const
   */
  CompressedReadSet& loadFromFile(const string& fileName);
};

#endif
