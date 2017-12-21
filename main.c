#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979
#define LU_tol 0.0001

// Structures

struct vector {
    size_t length;
    double *values;
};

struct matrix {
    size_t rows, columns;
    double **values;
};

// Constructors and destructors

struct vector *
vector_create(size_t length)
{
    struct vector *out = malloc(sizeof(*out));
    out->length = length;
    out->values = malloc(out->length * sizeof(out->values[0]));
    return out;
}

void
vector_destroy(struct vector *v)
{
    free(v->values);
    free(v);
}

struct matrix *
matrix_create(size_t rows, size_t cols)
{
    struct matrix *out = malloc(sizeof(*out));
    out->rows    = rows;
    out->columns = cols;
    out->values  = malloc(out->rows * sizeof(out->values[0]));
    for (size_t i = 0; i < out->rows; i++)
        out->values[i] = malloc(out->columns * sizeof(out->values[i][0]));
    return out;
}

void
matrix_destroy(struct matrix *m)
{
    for (size_t i = 0; i < m->rows; i++)
        free(m->values[i]);
    free(m->values);
    free(m);
}

// vector-only operations

double
vector_magnitude(struct vector *v)
{
    double out = 0;
    for (size_t i = 0; i < v->length; i++)
        out += (v->values[i] * v->values[i]);
    out = sqrt(out);
    return out;
}

struct vector *
vector_normal(struct vector *v)
{
    struct vector *out = vector_create(v->length);
    double mag = vector_magnitude(v);
    for (size_t i = 0; i < v->length; i++)
        out->values[i] = v->values[i] / mag;
    return out;
}

void
vector_print(struct vector *v)
{
    printf("[");
    for (size_t i = 0; i < v->length; i++)
        printf("%6.4f ", v->values[i]);
    printf("]\n");
}

struct vector *
vector_scalar_mult(struct vector *v, double mult)
{
    struct vector *out = vector_create(v->length);
    for (size_t i = 0; i < out->length; i++)
        out->values[i] = mult * v->values[i];
    return out;
}

struct vector *
vector_add(struct vector *v1, struct vector *v2)
{
    assert(v1->length == v2->length);
    struct vector *out = vector_create(v1->length);
    for (size_t i = 0; i < out->length; i++)
        out->values[i] = v1->values[i] + v2->values[i];
    return out;
}

struct vector *
vector_sub(struct vector *v1, struct vector *v2)
{
    assert(v1->length == v2->length);
    struct vector *temp = vector_scalar_mult(v2, -1.0);
    struct vector *out  = vector_add(v1, temp);
    vector_destroy(temp);
    return out;
}

double
vector_dot(struct vector *v1, struct vector *v2)
{
    assert(v1->length == v2->length);
    double out = 0;
    for (size_t i = 0; i < v1->length; i++)
        out += (v1->values[i] * v2->values[i]);
    return out;
}

double
vector_angle(struct vector *v1, struct vector *v2)
{
    assert(v1->length == v2->length);
    double dot = vector_dot(v1, v2);
    double m1  = vector_magnitude(v1);
    double m2  = vector_magnitude(v2);
    return acos(dot / (m1*m2));
}

struct vector *
vector_cross(struct vector *v1, struct vector *v2)
{
    assert(v1->length == v2->length);
    assert(v1->length == 3);
    struct vector *out = vector_create(v1->length);
    out->values[0] = (v1->values[1] * v2->values[2]) - (v1->values[2] * v2->values[1]);
    out->values[1] = (v1->values[2] * v2->values[0]) - (v1->values[0] * v2->values[2]);
    out->values[2] = (v1->values[0] * v2->values[1]) - (v1->values[1] * v2->values[0]);
    return out;
}

// matrix operations

void
matrix_print(struct matrix *m)
{
    printf("[");
    for (size_t i = 0; i < m->rows; i++) {
        if (i == 0)
            printf("[ ");
        else
            printf(" [ ");
        for (size_t j = 0; j < m->columns; j++)
            printf("%6.4f ", m->values[i][j]);
        printf("]");
        if (i != (m->rows - 1))
            printf("\n");
    }
    printf("]\n");
}

struct matrix *
matrix_clone(struct matrix *m)
{
    struct matrix *out = matrix_create(m->rows, m->columns);
    for (size_t i = 0; i < m->rows; i++)
        for (size_t j = 0; j < m->columns; j++)
            out->values[i][j] = m->values[i][j];
    return out;
}

struct matrix *
matrix_scalar_mul(struct matrix *m, double val)
{
    struct matrix *out = matrix_create(m->rows, m->columns);
    for (size_t i = 0; i < m->rows; i++)
        for (size_t j = 0; j < m->columns; j++)
            out->values[i][j] = val * m->values[i][j];
    return out;
}

struct matrix *
matrix_add(struct matrix *m1, struct matrix *m2)
{
    assert(m1->rows == m2->rows);
    assert(m1->columns == m2->columns);
    struct matrix *out = matrix_create(m1->rows, m1->columns);
    for (size_t i = 0; i < m1->rows; i++)
        for (size_t j = 0; j < m1->columns; j++)
            out->values[i][j] = m1->values[i][j] + m2->values[i][j];
    return out;
}


struct matrix *
matrix_sub(struct matrix *m1, struct matrix *m2)
{
    assert(m1->rows == m2->rows);
    assert(m1->columns == m2->columns);
    struct matrix *temp = matrix_scalar_mul(m2, -1.0);
    struct matrix *out  = matrix_add(m1, temp);
    matrix_destroy(temp);
    return out;
}

struct matrix *
matrix_mul(struct matrix *m1, struct matrix *m2)
{
    assert(m1->columns == m2->rows);
    struct matrix *out = matrix_create(m1->rows, m2->columns);
    for (size_t i = 0; i < out->rows; i++)
        for (size_t j = 0; j < out->columns; j++) {
            out->values[i][j] = 0;
            for (size_t k = 0; k < m1->columns; k++)
                out->values[i][j] += (m1->values[i][k] * m2->values[k][j]);
        }
    return out;
}

int
matrix_LU_decompose(struct matrix *m, int *P)
{
    assert(m->rows == m->columns);
    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i < m->rows; i++)
        P[i] = i;
    for (i = 0; i < m->rows; i++) {
        maxA = 0.0;
        imax = i;
        for (k = i; k < m->rows; k++)
            if ((absA = fabs(m->values[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }
        if (maxA < LU_tol)
            return 0;
        if (imax != i) {
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;
            
            ptr = m->values[i];
            m->values[i] = m->values[imax];
            m->values[imax] = ptr;

            P[m->rows]++;
        }

        for (j = i + 1; j < m->rows; j++) {
            m->values[j][i] /= m->values[i][i];

            for (k = i + 1; k < m->rows; k++)
                m->values[j][k] -= (m->values[j][i] * m->values[i][k]);
        }
    }
    return 1;
}

double
matrix_det(struct matrix *m)
{
    struct matrix *LU = matrix_clone(m);
    int *P            = malloc((1 + m->rows) * sizeof(P[0]));
    
    int test = matrix_LU_decompose(LU, P);
    
    if (!test) {
        printf("Decomposition failed.\n");
        return 0;
    }

    double det = LU->values[0][0];

    for (size_t i = 1; i < LU->rows; i++)
        det *= LU->values[i][i];

    if (((P[LU->rows] - LU->rows) % 2) != 0)
        det = -det;
    free(P);
    matrix_destroy(LU);
    
    return det;
}

struct matrix *
matrix_invert(struct matrix *m)
{
    struct matrix *LU = matrix_clone(m);
    int *P            = malloc((1 + m->rows) * sizeof(P[0]));

    int test = matrix_LU_decompose(LU, P);

    if (!test) {
        printf("Decomposition failed.\n");
        return NULL;
    }

    struct matrix *out = matrix_create(m->rows, m->columns);

    for (size_t j = 0; j < LU->rows; j++) {
        for (size_t i = 0; i < LU->rows; i++) {
            if (P[i] == j)
                out->values[i][j] = 1.0;
            else
                out->values[i][j] = 0.0;

            for (size_t k = 0; k < i; k++)
                out->values[i][j] -= LU->values[i][k] * out->values[k][j];
        }
        for (int i = (LU->rows - 1); i >= 0; i--) {
            for (int k = i + 1; k < LU->rows; k++)
                out->values[i][j] -= LU->values[i][k] * out->values[k][j];
            out->values[i][j] /= LU->values[i][i];
        }
    }
    
    free(P);
    matrix_destroy(LU);
    return out;
}

struct vector *
vector_matrix_mul(struct vector *v, struct matrix *m)
{
    assert(v->length == m->rows);
    struct vector *out = vector_create(m->columns);
    for (size_t i = 0; i < out->length; i++) {
        out->values[i] = 0;
        for (size_t j = 0; j < v->length; j++)
            out->values[i] += v->values[j] * m->values[j][i];
    }
    return out;
}

struct vector *
matrix_vector_mul(struct matrix *m, struct vector *v)
{
    assert(v->length == m->columns);
    struct vector *out = vector_create(m->rows);
    for (size_t i = 0; i < out->length; i++) {
        out->values[i] = 0;
        for (size_t j = 0; j < v->length; j++)
            out->values[i] += m->values[i][j] * v->values[j];
    }
    return out;
}

struct vector *
vector_LU_solve(struct matrix *m, struct vector *b)
{
    struct matrix *m_inv = matrix_invert(m);
    struct vector *out = matrix_vector_mul(m_inv, b);
    matrix_destroy(m_inv);
    return out;
}

// testbed

int main(void)
{
    struct vector *v1 = vector_create(3);
    struct vector *v2 = vector_create(3);
    v1->values[0] = 1;
    v1->values[1] = 0;
    v1->values[2] = 0;
    v2->values[0] = 0;
    v2->values[1] = 1;
    v2->values[2] = 0;
    struct vector *v3 = vector_cross(v1, v2);
    vector_destroy(v1);
    vector_destroy(v2);
    vector_destroy(v3);
    return 0;
}
