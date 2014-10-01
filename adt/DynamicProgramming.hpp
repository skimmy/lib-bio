#ifndef DYNAMIC_PROGRAMMING_H
#define DYNAMIC_PROGRAMMING_H

#include <cstddef>

typedef struct {// _matrix_point_2d {
  /**
   * Row index
   */
  size_t i;
  /**
   * Column index
   */
  size_t j;
} MatrixPoint2D;

//typedef MatrixPoint2D _matrix_point_2d;

/**
 * \brief This class is used to represent the general structure
 * of a dynamic programing model.
 *
 * This is an abstract class and maintains only information about
 * the dimension of the <em>Dynamic Programming</em> matrix. The
 * actual implementation of the matrix is left to the subclasses
 * to allow the maximum flexibility while obeying a common interface
 * to all those dynamic programming algorithms 
 *
 * The class defines several \c virtual methods that subclasses must
 * implement in order to provide the proper dymanic programming
 * algorithm.
 */
class DynamicProgramming {
protected:
  /**
   * \brief The number of rows of the dynamic programming matrix 
   */ 
  size_t n;
  /**
   * \brief The number of columns of the dynamic programming matrix 
   */ 
  size_t m;
public:
  // ---------------------------------------------------------
  //                       CONSTRUCTORS
  // ---------------------------------------------------------
  /**
   * \brief Creates an empty dynamic programming structure
   *
   * Default constructor simple create and empty dynamic 
   * programing strucutre (i.e. the number of rows and columns
   * for the DP matrix are both zero).
   *
   * \sa DynamicProgramming(size_t n, size_t m)
   */
  DynamicProgramming() { this->n = 0; this->m = 0; }
  /**
   * \brief Creates a dynamic programming structure with the
   * given dimensions
   *
   * \sa DynamicProgramming()
   */
  DynamicProgramming(size_t n, size_t m) { this->n = n; this->m = m; }

  // ---------------------------------------------------------
  //                        GET METHODS
  // ---------------------------------------------------------
  /**
   * \brief Returns the number of rows of the structure
   *
   * \return The number of rows
   *
   * \sa getColumnCount()
   */
  size_t getRowCount() { return n; }
  /**
   * \brief Returns the number of columns of the structure
   *
   * \return The number of columns
   *
   * \sa  getRowCount()
   */
  size_t getColumnCount() { return m; }

  // ---------------------------------------------------------
  //                      VIRTUAL METHODS
  // ---------------------------------------------------------
  /**
   * \brief Initializes the dynamic programming matrix
   *
   * The implementation of this method should take care of
   * initializing the dynamic programming matrix with the proper
   * values. All information other than the matrix dimensions
   * must be available at the subclass.
   *
   * \sa computeEntry(size_t i, size_t j)
   * \sa computeMatrix()
   */
  virtual void initMatrix() = 0;
  /**
   * \brief Computes an entry of the matrix.
   *
   * The implementation of this method should provide the proper
   * code to compute the entry \f$ (i,j) \f$ of the dynamic
   * programming matrix.
   *
   * \param i The row index of the entry
   * \param j The column index of the entry
   *
   * \sa initMatrix()
   * \sa computeMatrix()
   */
  virtual void computeEntry(size_t i, size_t j) = 0;
  /**
   * \brief Compute the dynamic programming matrix
   *
   * The implementation of this method should provide the proper
   * code to fill the whole dynamic programming matrix, it should
   * rely on other methods like initMatrix() and computeEntry(), 
   * or it should provide the whole filling code (depending on the
   * subclass design)
   *
   * \sa computeEntry(size_t i, size_t j)
   * \sa initMatrix()
   */
  virtual void computeMatrix() = 0;
  /**
   * \brief Prints the dynamic programming structure
   *
   * This virtual method is provided fo convinience and should
   * be used only for debug purposes.
   */
  virtual void printMatrix() = 0;
};

#endif
