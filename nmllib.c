//Nate Hall - Linear Algebra Library Tutorial Output

//This file contains the main matrix library, including the outcome of the tutorial
//All exceptions defined were pulled from the tutorial author's GitHub page
//Everything after a certain point was done by myself according to his tutorial

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "nmllib.h"
#include "nml_util.c"

#define DEFAULT_VALUE 0.0

#define CANNOT_ADD "Cannot add two matrices with different dimensions.\n"

#define CANNOT_SUBTRACT "Cannot subctract two matrices with different dimensions.\n"

#define CANNOT_MULITPLY \
  "Cannot multiply two matrices where \
  the number of columns of the first one \
  is different than the number of rows of the second one.\n" \

#define CANNOT_REMOVE_COLUMN "Cannot remove matrix column %d. The value should be less than %d.\n" 

#define CANNOT_REMOVE_ROW "Cannot remove matrix row %d. The value should be less than %d.\n" 

#define INVALID_ROWS \
  "Cannot create matrix with 0 number of rows. Aborting.\n" \

#define INVALID_COLS \
    "Cannot create matrix with 0 number of cols. Aborting.\n" \

#define CANNOT_TRACE \
    "Cannot calculate trace. Matrix needs to be square.\n" \

#define CANNOT_CROUT \
    "Cannot apply crout algorithm. Matrix needs to be square.\n" \

#define CANNOT_SWAP_ROWS \
     "Cannot swap rows (%d, %d) because the matrix number of rows is %d.\n" \

#define CANNOT_SWAP_COLUMNS \
      "Cannot swap columns (%d, %d) because the matrix number or columns is %d.\n" \

#define CANNOT_ROW_MULTIPLY \
      "Cannot multiply row (%d), maximum number of rows is %d.\n" \

#define CANNOT_COL_MULTIPLY "Cannot multiply col (%d), maximum number of columns is %d.\n" 
  
#define CANNOT_ADD_TO_ROW \
      "Cannot add %2.2f x (row=%d) to row=%d. Total number of rows is: %d.\n" \

#define CANNOT_LU_MATRIX_SQUARE \
      "Canot LU. Matrix (%d, %d) needs to be square.\n" \

#define CANNOT_LU_MATRIX_DEGENERATE \
      "Cannot LU. Matrix is degenerate or almost degenerate.\n" \

#define CANNOT_SOLVE_LIN_SYS_INVALID_B \
      "Cannot solve system. b[%d][%d] should have size b[%d][%d].\n" \

#define CANNOT_SET_DIAG \
      "Cannot set diag with value(=%2.2f). Matrix is not square.\n" \

#define CANNOT_CONCATENATE_H \
      "Cannot concatenate. Matrices have a different number of rows. Expected %d, found: %d.\n" \

#define CANNOT_CONCATENATE_V \
      "Cannot concatenate. Matrices have a different number of cols. Expected %d, found: %d.\n" \

#define CANNOT_GET_COLUMN \
      "Cannot get column (%d). The matrix has %d number of columns.\n" \

#define CANNOT_GET_ROW \
      "Cannot get row (%d). The matrix has %d number of rows.\n" \

#define INCONSITENT_ARRAY \
      "Cannot found element %d in the array (NULL). Expected a total of : %d elements.\n"  \

#define INCONSITENT_VARGS \
      "Cannot find element %d in the varargs. Expecteda total of : %d varargs.\n" \

#define CANNOT_REF_MATRIX_DEGENERATE \
      "Cannot compute REF. Matrix is degenerate or near degenerate.\n" \

#define CANNOT_OPEN_FILE "Cannot open file '%s'. Please check the path is correct and you have reading rights.\n"

#define INVALID_MATRIX_FILE \
      "Invalid matrix file: %s. Cannot read data.\n" \

#define VECTOR_J_DEGENERATE \
      "Vector on colum %d is generate or near degenerate. Cannot proceed further.\n" \

#define CANNOT_QR_NON_SQUARE \
      "We cannot QA non-square matrix[%d, %d].\n" \

#define CANNOT_COLUMN_L2NORM \
      "Cannot get column (%d). The matrix has %d numbers of columns.\n" \

#define CANNOT_VECT_DOT_DIMENSIONS \
      "The two vectors have different dimensions: %d and %d.\n" \
       

//Everything after this section is from me following the tutorial, with my comments
//*****************************************************************************
//Matrix creation section

//struct that represents our matrix. declared in nmllib.h
/*
typedef struct nml_mat_s {
  unsigned int num_rows;
  unsigned int num_cols;
  double **data;
  int is_square;
} nml_mat;
*/

nml_mat *nml_mat_new(unsigned int num_rows, unsigned int num_cols){
    //input validation for zero rows
    if (num_rows == 0) {
        NML_ERROR(INVALID_ROWS);
        return NULL;
    }
    if (num_cols == 0) {
        NML_ERROR(INVALID_COLS);
        return NULL;
    }
    //here we allocate memory for a basic nml_mat struct
    nml_mat *m = calloc(1, sizeof(*m));
    NP_CHECK(m);
    m->num_rows = num_rows;
    m->num_cols = num_cols;
    m->is_square = (num_rows == num_cols) ? 1:0;
    //now we begin allocating memory for the matrix itself, first by row, then column
    m->data = calloc(m->num_rows, sizeof(*m->data));
    NP_CHECK(m->data);
    int i;
    for(i=0; i < m->num_rows; i++) {
        m->data[i] = calloc(m->num_cols, sizeof(**m->data));
        NP_CHECK(m->data[i]);
    }
    return m;
}


void nml_mat_free(nml_mat *matrix) {
    //deallocate memory row by row, then the list of pointers, then the struct itself
    int i;
    for(i = 0; i < matrix->num_rows; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);
}


nml_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, double min, double max) {
    nml_mat *r = nml_mat_new(num_rows, num_cols);
    int i, j;
    for(i=0; i < num_rows; i++) {
        for(j=0; j < num_cols; j++) {
            r->data[i][j] = nml_rand_interval(min,max);
            //generate a random number for each matrix entry
            //the nml_rand_interval method is defined in nml_util.c
        }
    }
    return r;
}

//constructs a square matrix
nml_mat *nml_mat_sqr(unsigned int size) {
    return nml_mat_new(size,size);
}

//constructs a random square matrix
nml_mat *nml_mat_sqr_rnd(unsigned int size, double min, double max) {
    return nml_mat_rnd(size, size, min, max);
}

//constructs the identity matrix with *size* rows and *size* columns
nml_mat *nml_mat_eye(unsigned int size) {
    nml_mat *r = nml_mat_new(size,size);
    int i;
    for(i=0; i < r->num_rows; i++) {
        r->data[i][i] = 1.0;
    }
    return r;
}

//copies the data from one matrix to another
nml_mat *nml_mat_cp(nml_mat *m) {
    nml_mat *r = nml_mat_new(m->num_rows,m->num_cols);
    int i,j;
    for(i=0; i < m->num_rows; i++) {
        for(j=0; j < m->num_cols; j++) {
            r->data[i][j] = m->data[i][j];
        }
    }
    return r;
}

//Reads input from a file (with certain formatting) and creates a matrix
//Can also take keyboard input by calling nml_mat_fromfilef(stdin);
nml_mat *nml_mat_fromfilef(FILE *f) {
    int i, j;
    unsigned int num_rows = 0, num_cols = 0;
    fscanf(f, "%d", &num_rows);
    fscanf(f, "%d", &num_cols);
    nml_mat *r = nml_mat_new(num_rows, num_cols);
    for(i=0; i < num_rows; i++) {
        for(j=0; j< num_cols; j++) {
            fscanf(f, "%lf\t", &r->data[i][j]);
        }
    }
    return r;
}

//*****************************************************************************
//Matrix method section

//Checks to see if matrices have equal dimension
int nml_mat_eqdim(nml_mat *m1, nml_mat *m2) {
    return (m1->num_cols == m2->num_cols) && (m1->num_rows == m2->num_rows);
}

//Checks for matrix equality given a certain tolerance
//For exact equality use tolerance = 0.0
int nml_mat_eq(nml_mat *m1, nml_mat *m2, double tolerance) {
    //checks dimension equality first
    if (!nml_mat_eqdim(m1, m2)) {
        return 0;
    }
    int i,j;
    for(i = 0; i < m1->num_rows; i++) {
        for(j = 0; j < m1->num_cols; j++) {
            //then for each element in m1, check if the corresponding element in m2
            //is within *tolerance* of m1[i][j]
            if(fabs(m1->data[i][j] - m2->data[i][j]) > tolerance) {
                return 0;
            }
        } 
    }
    return 1;
}

//Prints the matrix on the stdout (with a custom formatting for elements)
void nml_mat_printf(nml_mat *matrix, const char *d_fmt) {
    int i,j;
    fprintf(stdout, "\n");
    for(i = 0; i < matrix->num_rows; ++i) {
        for(j = 0; j < matrix->num_cols; ++j) {
            fprintf(stdout, d_fmt, matrix->data[i][j]);
        }
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}

//Special print case with "default" formatting
void nml_mat_print(nml_mat *matrix) {
  nml_mat_printf(matrix, "%lf\t\t");
}

//Method that retrieves a column from a matrix- will be useful in other methods
nml_mat *nml_mat_col_get(nml_mat *m, unsigned int col) {
    if (col >= m->num_cols) {
        NML_FERROR(CANNOT_GET_COLUMN, col, m->num_cols);
        return NULL;
    }
    //Create a new 1-column matrix to hold the retrieved elements
    //then fill it with the elements from the given matrix at the given column
    nml_mat *r = nml_mat_new(m->num_rows, 1);
    int j;
    for(j = 0; j < r->num_rows; j++) {
        r->data[j][0] = m->data[j][col];
    }
    return r;
}

//Similar method that retrieves a row from a matrix
//Simpler to implement than columns since we use row-first indexing
nml_mat *nml_mat_row_get(nml_mat *m, unsigned int row) {
    if (row >= m->num_rows) {
        NML_FERROR(CANNOT_GET_ROW, row, m->num_rows);
        return NULL;
    }
    nml_mat *r = nml_mat_new(1, m->num_cols);
    memcpy(r->data[0], m->data[row], m->num_cols * sizeof(*r->data[0]));
    return r;
} 

//helper method to set all the elements of a matrix to a given value
void nml_mat_all_set(nml_mat *matrix, double value) {
    int i, j;
    for(i = 0; i < matrix->num_rows; i++) {
        for(j = 0; j < matrix->num_cols; j++) {
            matrix->data[i][j] = value;
        }
    }
}

//helper method to set all the elements of a matrix's diagonal to a given value
int nml_mat_diag_set(nml_mat *m, double value) {
    //note - in implementing SVD it may be useful to be able to
    //generalize this method to matrices of arbitrary size
    if (!m->is_square) {
        NML_FERROR(CANNOT_SET_DIAG, value);
        return 0;
    }
    int i;
    for(i = 0; i < m->num_rows; i++) {
        m->data[i][i] = value;
    }
    return 1;
} 

//Method that multiplies a row of a given matrix by a given scalar
//In-place implementation
int nml_mat_row_mult_r(nml_mat *m, unsigned int row, double num) {
    if (row>= m->num_rows) {
        NML_FERROR(CANNOT_ROW_MULTIPLY, row, m->num_rows);
        return 0;
    }
    int i;
    for(i=0; i < m->num_cols; i++) {
        m->data[row][i] *= num;
    }
    return 1;
}

//Method that multiplies a row of a given matrix by a given scalar
//Retrieves a copy of the matrix with the given row multiplication
nml_mat *nml_mat_row_mult(nml_mat *m, unsigned int row, double num) {
    nml_mat *r = nml_mat_cp(m);
    if (!nml_mat_row_mult_r(r, row, num)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
}

//Method that multiplies a column of a given matrix by a given scalar
//In-place implementation
int nml_mat_col_mult_r(nml_mat *m, unsigned int col, double num) {
    if (col>=m->num_cols) {
        NML_FERROR(CANNOT_COL_MULTIPLY, col, m->num_cols);
        return 0;
    }
    int i;
    for(i = 0; i < m->num_rows; i++) {
        m->data[i][col] *= num;
    }
    return 1;
} 

//Method that multiplies a column of a given matrix by a given scalar
//Retrieves a copy of the matrix with the given column multiplication
nml_mat *nml_mat_col_mult(nml_mat *m, unsigned int col, double num) {
    nml_mat *r = nml_mat_cp(m);
    if (!nml_mat_col_mult_r(r, col, num)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
}

//Two methods to add a constant multiple of a row to another row, 
//an elementary row operation.
//Similarly to row and column multiplication, we will have
//both an in-place and a copy retrieval implementation.

//Adds a multiple of the row with index *row* to the row with index *where*
int nml_mat_row_addrow_r(nml_mat *m, unsigned int where, unsigned int row, double multiplier) {
    if (where >= m->num_rows || row >= m->num_rows) {
        NML_FERROR(CANNOT_ADD_TO_ROW, multiplier, row, where, m->num_rows);
        return 0;
    }
    int i = 0;
    for(i = 0; i < m->num_cols; i++) {
        m->data[where][i] += multiplier * m->data[row][i];
    }
    return 1;
} 

//Returns a copy of the given matrix with a multiple of
//the row with index *row* added to the row with index *where*
nml_mat *nml_mat_row_addrow(nml_mat *m, unsigned int where, unsigned int row, double multiplier) {
    nml_mat *r = nml_mat_cp(m);
    if (!nml_mat_row_addrow_r(m, where, row, multiplier)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
}

//Multiplies an entire matrix by a given scalar, in-place
int nml_mat_smult_r(nml_mat *m, double num) {
    int i, j;
    for(i = 0; i < m->num_rows; i++) {
        for(j = 0; j < m->num_cols; j++) {
            m->data[i][j] *= num;
        }
    }
    return 1;
} 

//Returns a copy of the given matrix multiplied by a given scalar
nml_mat *nml_mat_smult(nml_mat *m, double num) {
    nml_mat *r = nml_mat_cp(m);
    nml_mat_smult_r(r, num);
    return r;
}

//Method that returns a copy of a given matrix with a specified column removed
nml_mat *nml_mat_col_rem(nml_mat *m, unsigned int column) {
    if(column >= m->num_cols) {
        NML_FERROR(CANNOT_REMOVE_COLUMN, column, m->num_cols);
        return NULL;
    }
    //create new matrix with different dimensions
    nml_mat *r = nml_mat_new(m->num_rows, m->num_cols-1);
    //then copy over the relevant data, keeping different indices for
    //the columns of r and the columns of m
    int i, j, k;
    for(i = 0; i < m->num_rows; i++) {
        //line below originally set k=0, but this may cause a bug if k
        //is immediately incremented, causing the first column to be missed.
        //testing required.
        for(j = 0, k=-1; j < m->num_cols; j++) {
            if (column!=j) {
                r->data[i][k++] = m->data[i][j];
            }
        }
    }
  return r;
}

//Method that returns a copy of a given matrix with a specified row removed
nml_mat *nml_mat_row_rem(nml_mat *m, unsigned int row) {
    if (row >= m->num_rows) {
        NML_FERROR(CANNOT_REMOVE_ROW, row, m->num_rows);
        return NULL;
    }
    //similar to last method, with less danger of skipping a necessary row
    nml_mat *r = nml_mat_new(m->num_rows-1, m->num_cols);
    int i, j, k;
    for(i = 0, k = 0; i < m->num_rows; i++) {
        if (row!=i) {
            for(j = 0; j < m->num_cols; j++) {
                r->data[k][j] = m->data[i][j];
            }
            k++;
        }
    }
  return r;
}

//Method that swaps two specified rows of a matrix in-place
int nml_mat_row_swap_r(nml_mat *m, unsigned int row1, unsigned int row2) {
    if (row1 >= m->num_rows || row2 >= m->num_rows) {
        NML_FERROR(CANNOT_SWAP_ROWS, row1, row2, m->num_rows);
        return 0;
    }
    //note that the pointer implementation makes this very simple
    double *tmp = m->data[row2];
    m->data[row2] = m->data[row1];
    m->data[row1] = tmp;
    return 1;
} 

//Method that returns a copy of a given matrix with two specified rows swapped
nml_mat *nml_mat_row_swap(nml_mat *m, unsigned int row1, unsigned int row2) {
    nml_mat *r = nml_mat_cp(m);
    if (!nml_mat_row_swap_r(r, row1, row2)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
} 

//Method that swaps two specified columns of a given matrix in-place
int nml_mat_col_swap_r(nml_mat *m, unsigned int col1, unsigned int col2) {
    if (col1 >= m->num_cols || col2 >= m->num_cols) {
        NML_FERROR(CANNOT_SWAP_ROWS, col1, col2, m->num_cols);
        return 0;
    }
    //Since matrices are implemented using row-first indexing,
    //this is not as simple or efficient as swapping rows
    double tmp;
    int j;
    for(j = 0; j < m->num_rows; j++) {
        tmp = m->data[j][col1];
        m->data[j][col1] = m->data[j][col2];
        m->data[j][col2] = tmp;
    }
    return 1;
}

//Method that returns a copy of the given matrix with two specified columns swapped
nml_mat *nml_mat_col_swap(nml_mat *m, unsigned int col1, unsigned int col2) {
    nml_mat *r = nml_mat_cp(m);
    if (!nml_mat_col_swap_r(r, col1, col2)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
} 

//Original tutorial here implements methods to horizontally and vertically
//concatenate an array of matrices.
//I will skip this for now, but in the future may implement a special case of
//these operations to reconstruct a matrix out of block matrices,
//as well as break a matrix into block matrices,
//as these operations may be useful in implementing the Implicit QR Step algorithm.

//*****************************************************************************
//Matrix operation section

//Method that adds two matrices together, directly editing the entries of m1
int nml_mat_add_r(nml_mat *m1, nml_mat *m2) {
    if (!nml_mat_eqdim(m1, m2)) {
        NML_ERROR(CANNOT_ADD);
        return 0;
    }
    int i, j;
    for(i = 0; i < m1->num_rows; i++) {
        for(j = 0; j < m1->num_rows; j++) {
            m1->data[i][j] += m2->data[i][j];
            }
    }
    return 1;
}

//Method that returns a new matrix that is the sum of given matrices m1 and m2
nml_mat *nml_mat_add(nml_mat *m1, nml_mat *m2) {
    nml_mat *r = nml_mat_cp(m1);
    if (!nml_mat_add_r(r, m2)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
}

//Next two matrix operations implement subtracting one matrix from another
//Including these since they are in the original tutorial,
//but in the future I may combine these both into a more general method
//that adds together scalar multiples of given matrices

//Subtracts m2 from m1, directly editing the entries of m1
int nml_mat_sub_r(nml_mat *m1, nml_mat *m2) {
    if (!nml_mat_eqdim(m1, m2)) {
        NML_ERROR(CANNOT_SUBTRACT);
        return 0;
    }
    int i, j;
    for(i = 0; i < m1->num_rows; i++) {
        for(j = 0; j < m1->num_cols; j++) {
        m1->data[i][j] -= m2->data[i][j];
        }
    }
    return 1;
} 

//Returns a new matrix equal to m1 - m2
nml_mat *nml_mat_sub(nml_mat *m1, nml_mat *m2) {
    nml_mat *r = nml_mat_cp(m2);
    if (!nml_mat_sub_r(r, m2)) {
        nml_mat_free(r);
        return NULL;
    }
    return r;
}

//Method that implements naive matrix multiplication
//Returns a new matrix equal to (m1)(m2)
nml_mat *nml_mat_dot(nml_mat *m1, nml_mat *m2) {
    if (!(m1->num_cols == m2->num_rows)) {
        NML_ERROR(CANNOT_MULITPLY);
        return NULL;
    }
    int i, j, k;
    nml_mat *r = nml_mat_new(m1->num_rows, m2->num_cols);
    for(i = 0; i < r->num_rows; i++) {
        for(j = 0; j < r->num_cols; j++) {
            for(k = 0; k < m1->num_cols; k++) {
                r->data[i][j] += m1->data[i][k] * m2->data[k][j];
            }
        }
    }
    return r;
}

//Helper method for implementing row reduction
//Returns the index of a pivot, i.e. the first nonzero element in the *col*th column,
//below the *row*th row. If no pivot is found, returns -1
int _nml_mat_pivotidx(nml_mat *m, unsigned int col, unsigned int row) {
    // No validations are made, as this is a helper method
    int i;
    for(i = row; i < m->num_rows; i++) {
        //checks entry against a library constant NML_MIN_COEF
        //instead of zero to prevent overflow
        if (fabs(m->data[i][col]) > NML_MIN_COEF) {
            return i;
        }
    }
    return -1;
}

//Returns the matrix in Row Echelon Form, 
//which is found by using row reduction/Gaussian elimination
nml_mat *nml_mat_ref(nml_mat *m) {
    nml_mat *r = nml_mat_cp(m);
    int i, j, k, pivot;
    j = 0, i = 0;
    //We iterate until we exhaust the columns and the rows
    while(j < r->num_cols && i < r->num_cols) { //should this say i < r->num_rows?
        // Find the pivot
        pivot = _nml_mat_pivotidx(r, j, i);
        if (pivot<0) {
            //All elements below row i on the column are zeros so we move to the next column
            j++;
            continue;
        }
        //We move the pivot to the first row that doesn't a pivot in place already
        if (pivot!=i) {
            nml_mat_row_swap_r(r, i, pivot);
        }
        //Multiply each element in the pivot row by the reciprocal of the pivot
        //i.e. normalize the pivot to create a leading 1
        nml_mat_row_mult_r(r, i, 1/r->data[i][j]);
        //Add multiples of the pivot so every element below the pivot in column j is equal to 0
        for(k = i+1; k < r->num_rows; k++) {
            //Check against NML_MIN_COEF similarly to the pivotidx method
            //We consider anything below NML_MIN_COEF to be 0
            //Note to self: might we want to cement this somehow to make the method airtight?
            if (fabs(r->data[k][j]) > NML_MIN_COEF) {
                nml_mat_row_addrow_r(r, k, i, -(r->data[k][j]));
            } 
        }
        i++;
        j++;
    }
    return r;
} 

//Helper method for implementing Gauss-Jordan Elimination
//More numerically stable than our previous pivot helper method
//Instead of returning the index of the first nonzero element of the given column
//under the row, returns the maximum element in the same range.
//If no pivot is found, return -1
int _nml_mat_pivotmaxidx(nml_mat *m, unsigned int col, unsigned int row) {
    int i, maxi;
    double micol;
    double max = fabs(m->data[row][col]);
    maxi = row;
    for(i = row; i < m->num_rows; i++) {
        micol = fabs(m->data[i][col]);
        if (micol>max) {
            max = micol;
            maxi = i;
        }
    }
    return (max < NML_MIN_COEF) ? -1 : maxi;
} 

//Method that returns the given matrix in Reduced Row Echelon Form
nml_mat *nml_mat_rref(nml_mat *m) {
    nml_mat* r = nml_mat_cp(m);
    int i,j,k,pivot;
    i = 0;
    j = 0;
    while(j < r->num_cols && i < r->num_rows) {
        //Find the pivot, as detailed in the last method
        pivot = _nml_mat_pivotmaxidx(r, j, i);
        if (pivot<0) {
            // No pivot, we change columns
            j++;
            continue;
        }
        //Put the pivot row into the desired position
        if (pivot!=i) {
            nml_mat_row_swap_r(r, i, pivot);
        }
        //Normalize the leading coeff, create 1 in the pivot position
        nml_mat_row_mult_r(r, i, 1/r->data[i][j]);
        // We put zeros on the column with the pivot
        for(k = 0; k < r->num_rows; k++) {
            if (!(k==i)) {
                nml_mat_row_addrow_r(r, k, i, -(r->data[k][j]));
            }
        }
        i++;
        j++;
    }
    return r;
}

//*****************************************************************************
//LU(P) Decomposition

//The LUP struct defined in nmllib.h:
/*
typedef struct nml_mat_lup_s {
    nml_mat *L;
    nml_mat *U;
    nml_mat *P;
    unsigned int num_permutations;
} nml_mat_lup;
*/


//Here we create constructor- and destructor-like methods

nml_mat_lup *nml_mat_lup_new(nml_mat *L, nml_mat *U, nml_mat *P, unsigned int num_permutations) {
    nml_mat_lup *r = malloc(sizeof(*r));
    NP_CHECK(r);
    r->L = L;
    r->U = U;
    r->P = P;
    r->num_permutations = num_permutations;
    return r;
}

void nml_mat_lup_free(nml_mat_lup* lu) {
    nml_mat_free(lu->P);
    nml_mat_free(lu->L);
    nml_mat_free(lu->U);
    free(lu);
} 

// Finds the maxid on the column (starting from k -> num_rows)
// This method is used for pivoting in LUP decomposition
int _nml_mat_absmaxr(nml_mat *m, unsigned int k) {
    int i;
    double max = m->data[k][k];
    int maxIdx = k;
    for(i = k+1; i < m->num_rows; i++) {
        if (fabs(m->data[i][k]) > max) {
            max = fabs(m->data[i][k]);
            maxIdx = i;
        }
    }
    return maxIdx;
}

//Method that computes the LU(P) factorization and returns an nml_mat_lup struct
nml_mat_lup *nml_mat_lup_solve(nml_mat *m) {
    if (!m->is_square) {
        NML_FERROR(CANNOT_LU_MATRIX_SQUARE, m->num_rows, m-> num_cols);
        return NULL;
    }
    nml_mat *L = nml_mat_new(m->num_rows, m->num_rows);
    nml_mat *U = nml_mat_cp(m);
    nml_mat *P = nml_mat_eye(m->num_rows);

    int j,i, pivot;
    unsigned int num_permutations = 0;
    double mult;

    for(j = 0; j < U->num_cols; j++) {
        // Retrieves the row with the biggest element for column (j)
        pivot = _nml_mat_absmaxr(U, j);
        if (fabs(U->data[pivot][j]) < NML_MIN_COEF) {
            NML_ERROR(CANNOT_LU_MATRIX_DEGENERATE);
            return NULL;
        }
        if (pivot!=j) {
            //Pivots LU and P accordingly
            nml_mat_row_swap_r(U, j, pivot);
            nml_mat_row_swap_r(L, j, pivot);
            nml_mat_row_swap_r(P, j, pivot);
            //Keeps the number of permutations to easily calculate the determinant
            num_permutations++;
        }
        for(i = j+1; i < U->num_rows; i++) {
            mult = U->data[i][j] / U->data[j][j];
            //Building the U upper rows
            nml_mat_row_addrow_r(U, i, j, -mult);
            //Store the multiplier in L
            L->data[i][j] = mult;
        }
    }
    nml_mat_diag_set(L, 1.0);
    return nml_mat_lup_new(L, U, P, num_permutations);
} 

//The LU factorization can be used to solve systems of linear equations

//Forward substitution algorithm
//Solves the linear system L * x = b
//
//L is lower triangular matrix of size NxN
//B is column matrix of size Nx1
//x is the solution column matrix of size Nx1
//
//The algorithm still attempts to solve the system even if the matrix is
//not lower triangular- consider adding input validation
//
//Also note that the system is not solvable if any of the diagonal elements are 0
nml_mat *nml_ls_solvefwd(nml_mat *L, nml_mat *b) {
    nml_mat* x = nml_mat_new(L->num_cols, 1);
    int i,j;
    double tmp;
    for(i = 0; i < L->num_cols; i++) {
        tmp = b->data[i][0];
        for(j = 0; j < i ; j++) {
            tmp -= L->data[i][j] * x->data[j][0];
        }
        x->data[i][0] = tmp / L->data[i][i];
    }
    return x;
}

//Back substition algorithm
//Solves the linear system U *x = b
//
//U is an upper triangular matrix of size NxN
//B is a column matrix of size Nx1
//x is the solution column matrix of size Nx1
//
//Similarly to solvefwd, the algorithm still attempts to solve the system 
//even if the matrix is not lower triangular- consider adding input validation
//
//Also note that the system is not solvable if any of the diagonal elements are 0
nml_mat *nml_ls_solvebck(nml_mat *U, nml_mat *b) {
    nml_mat *x = nml_mat_new(U->num_cols, 1);
    int i = U->num_cols, j;
    double tmp;
    while(i-->0) {
        tmp = b->data[i][0];
        for(j = i; j < U->num_cols; j++) {
            tmp -= U->data[i][j] * x->data[j][0];
        }
        x->data[i][0] = tmp / U->data[i][i];
    }
    return x;
} 

//Uses the solvefwd and solvebck methods to solve an arbitrary (square) linear system
//Ax = b, where A's LU(P) factorization has already been computed
nml_mat *nml_ls_solve(nml_mat_lup *lu, nml_mat* b) {
    if (lu->U->num_rows != b->num_rows || b->num_cols != 1) {
        NML_FERROR(CANNOT_SOLVE_LIN_SYS_INVALID_B,
            b->num_rows,
            b->num_cols,
            lu->U->num_rows,
        1);
        return NULL;
    }
    nml_mat *Pb = nml_mat_dot(lu->P, b);

    //Solve L*y = P*b using forward substition
    nml_mat *y = nml_ls_solvefwd(lu->L, Pb);

    //Solve U*x=y using back substitution
    nml_mat *x = nml_ls_solvebck(lu->U, y);

    nml_mat_free(y);
    nml_mat_free(Pb);
    return x;
} 

//Calculates the inverse of an n by n matrix
//by solving n linear systems of n variables.
nml_mat *nml_mat_inv(nml_mat_lup *lup) {
    unsigned n = lup->L->num_cols;
    nml_mat *r = nml_mat_sqr(n);
    nml_mat *I = nml_mat_eye(lup->U->num_rows);
    nml_mat *invx;
    nml_mat *Ix;
    int i,j;
    for(j =0; j < n; j++) {
        Ix = nml_mat_col_get(I, j);
        invx = nml_ls_solve(lup, Ix);
        for(i = 0; i < invx->num_rows; i++) {
            r->data[i][j] = invx->data[i][0];
        }
        nml_mat_free(invx);
        nml_mat_free(Ix);
    }
    nml_mat_free(I);
    return r;
} 

//Calculating the determinant of a matrix is normally an O(n!) operation.
//However, after the LU(P) factorisation the determinant can be easily calculated
//by multiplying the main diagonal of matrix U with the sign.
//the sign is -1 if the number of permutations is odd
//the sign is +1 if the number of permutations is even
double nml_mat_det(nml_mat_lup* lup) {
    int k;
    int sign = (lup->num_permutations%2==0) ? 1 : -1;
    nml_mat *U = lup->U;
    double product = 1.0;
    for(k = 0; k < U->num_rows; k++) {
        product *= U->data[k][k];
    }
    return product * sign;
}

//*****************************************************************************
//QR Decomposition

//The QR struct defined in nmllib.h:
/*
typedef struct nml_mat_qr_s {
    nml_mat *Q;
    nml_mat *R;
} nml_mat_qr;
*/

//Some constructor- and destructor- like methods:

nml_mat_qr *nml_mat_qr_new() {
  nml_mat_qr *qr = malloc(sizeof(*qr));
  NP_CHECK(qr);
  return qr;
}

void nml_mat_qr_free(nml_mat_qr *qr) {
  nml_mat_free(qr->Q);
  nml_mat_free(qr->R);
  free(qr);
}


//Useful for QR decomposition
//Returns the dot product of two columns of two matrices
//May scrap later and create a slightly more general dot product method
double nml_vect_dot(nml_mat *m1, unsigned int m1col, nml_mat *m2, unsigned m2col) {
    if (m1->num_rows!=m2->num_rows) {
        NML_FERROR(CANNOT_VECT_DOT_DIMENSIONS, m1->num_rows, m2->num_rows);
    }
    if (m1col >= m1->num_cols) {
        NML_FERROR(CANNOT_GET_COLUMN, m1col, m1->num_cols);
    }
    if (m2col >= m2->num_cols) {
        NML_FERROR(CANNOT_GET_COLUMN, m2col, m2->num_cols);
    }
    int i;
    double dot = 0.0;
    for(i = 0; i < m1->num_rows; i++) {
        dot += m1->data[i][m1col] * m2->data[i][m2col];
    }
    return dot;
} 

//Calculates the l2 norm for a column in a given matrix
double nml_mat_col_l2norm(nml_mat *m, unsigned int col) {
    if (col >= m->num_cols) {
        NML_FERROR(CANNOT_COLUMN_L2NORM, col, m->num_cols);
    }
    double doublesum = 0.0;
    int i;
    for(i = 0; i < m->num_rows; i++) {
        doublesum += (m->data[i][col]*m->data[i][col]);
    }
    return sqrt(doublesum);
}

//Calculates the l2norm for each column of a given matrix
//Returns results in a 1-row matrix
nml_mat *nml_mat_l2norm(nml_mat *m) {
    int i, j;
    nml_mat *r = nml_mat_new(1, m->num_cols);
    double square_sum;
    for(j = 0; j < m->num_cols; j++) {
        square_sum = 0.0;
        for(i = 0; i < m->num_rows; i++) {
            square_sum+=m->data[i][j]*m->data[i][j];
        }
        r->data[0][j] = sqrt(square_sum);
    }
    return r;
} 

//Helper method for calculating the QR decomposition
//Normalizes the columns of a matrix, important for constructing orthogonal matrices
int nml_mat_normalize_r(nml_mat *m) {
  nml_mat *l2norms = nml_mat_l2norm(m);
  int j;
  for(j = 0; j < m->num_cols; j++) {
    if (l2norms->data[0][j] < NML_MIN_COEF) {
      NML_FERROR(VECTOR_J_DEGENERATE, j);
      nml_mat_free(l2norms);
      return 0;
    }
    nml_mat_col_mult_r(m, j, 1/l2norms->data[0][j]);
  }
  nml_mat_free(l2norms);
  return 1;
}

//Method that calculates the QR Decomposition of a given matrix
//according to the Gram-Schmidt process for finding an orthonormal basis
nml_mat_qr *nml_mat_qr_solve(nml_mat *m) {

    nml_mat_qr *qr = nml_mat_qr_new();
    nml_mat *Q = nml_mat_cp(m);
    nml_mat *R = nml_mat_new(m->num_rows, m->num_cols);

    int j, k;
    double l2norm;
    double rkj;
    nml_mat *aj;
    nml_mat *qk;
    for(j=0; j < m->num_cols; j++) {    
        rkj = 0.0;
        aj = nml_mat_col_get(m, j);
        for(k = 0; k < j; k++) {
            rkj = nml_vect_dot(m, j, Q, k);
            R->data[k][j] = rkj;
            qk = nml_mat_col_get(Q, k);
            nml_mat_col_mult_r(qk, 0, rkj);
            nml_mat_sub_r(aj, qk);
            nml_mat_free(qk);
        }
        for(k = 0; k < Q->num_rows; k++) {
            Q->data[k][j] = aj->data[k][0];
        }
        l2norm = nml_mat_col_l2norm(Q, j);
        nml_mat_col_mult_r(Q, j, 1/l2norm);
        R->data[j][j] = l2norm;
        nml_mat_free(aj);
    }
    qr->Q = Q;
    qr->R = R;
    return qr;
} 

//main method, so the compiler doesn't yell at me
int main() {

}