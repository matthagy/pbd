/*********************************************************************          
 *               Simple 2D and 3D Vector Library
 *--------------------------------------------------------------------
 * Copyright (C) 2009 Matthew Hagy <hagy@case.edu>
 *--------------------------------------------------------------------
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef Vec_SQRT_FUNC
#  define Vec_SQRT_FUNC sqrt
#endif

#ifndef Vec_ACOS_FUNC
#  define Vec_ACOS_FUNC acos 
#endif

#define Vec2_COMPONENTS(TYPE) TYPE x, y

typedef struct {
        Vec2_COMPONENTS(double);
} Vec2D;

typedef struct {
        Vec2_COMPONENTS(float);
} Vec2F;

typedef struct {
        Vec2_COMPONENTS(long);
} Vec2L;

typedef struct {
        Vec2_COMPONENTS(unsigned long);
} Vec2UL;

typedef struct {
        Vec2_COMPONENTS(int);
} Vec2I;

typedef struct {
        Vec2_COMPONENTS(unsigned int);
} Vec2UI;

typedef struct {
        Vec2_COMPONENTS(short);
} Vec2S;

typedef struct {
        Vec2_COMPONENTS(unsigned short);
} Vec2US;

typedef struct {
        Vec2_COMPONENTS(char);
} Vec2C;

typedef struct {
        Vec2_COMPONENTS(unsigned char);
} Vec2UC;


#define Vec3_COMPONENTS(TYPE) TYPE x,y,z

typedef struct {
        Vec3_COMPONENTS(double);
} Vec3D;

typedef struct {
        Vec3_COMPONENTS(float);
} Vec3F;

typedef struct {
        Vec3_COMPONENTS(long);
} Vec3L;

typedef struct {
        Vec3_COMPONENTS(unsigned long);
} Vec3UL;

typedef struct {
        Vec3_COMPONENTS(int);
} Vec3I;

typedef struct {
        Vec3_COMPONENTS(unsigned int);
} Vec3UI;

typedef struct {
        Vec3_COMPONENTS(short);
} Vec3S;

typedef struct {
        Vec3_COMPONENTS(unsigned short);
} Vec3US;

typedef struct {
        Vec3_COMPONENTS(char);
} Vec3C;

typedef struct {
        Vec3_COMPONENTS(unsigned char);
} Vec3UC;



#define Vec2_FRMT(FRMT) "<" FRMT " " FRMT ">"
#define Vec3_FRMT(FRMT) "<" FRMT " " FRMT " " FRMT ">"

#define Vec2_IFRMT(FRMT) "<" FRMT "." FRMT ">"
#define Vec3_IFRMT(FRMT) "<" FRMT "." FRMT "." FRMT ">"


#define Vec2D_FRMT Vec2_FRMT("%f")
#define Vec3D_FRMT Vec3_FRMT("%f")
#define Vec2F_FRMT Vec2_FRMT("%f")
#define Vec3F_FRMT Vec3_FRMT("%f")
#define Vec2E_FRMT Vec2_FRMT("%e")
#define Vec3E_FRMT Vec3_FRMT("%e")
#define Vec2G_FRMT Vec2_FRMT("%g")
#define Vec3G_FRMT Vec3_FRMT("%g")
#define Vec2I_FRMT Vec2_IFRMT("%i")
#define Vec3I_FRMT Vec3_IFRMT("%i")

#define Vec2_ARGS(VEC) (VEC).x, (VEC).y
#define Vec3_ARGS(VEC) (VEC).x, (VEC).y, (VEC).z
#define Vec3_ARGS_SCALED(SCALE, VEC) (SCALE)*(VEC).x, (SCALE)*(VEC).y, (SCALE)*(VEC).z

#define Vec2_SET(DEST, V1, V2)                                                 \
        (DEST).x = (V1),                                                       \
        (DEST).y = (V2)                                                        \

#define Vec3_SET(DEST, V1, V2, V3)                                             \
        (DEST).x = (V1),                                                       \
        (DEST).y = (V2),                                                       \
        (DEST).z = (V3)                                                        \

#define Vec3_SETALL(DEST, OP) do {              \
        __typeof__(OP) _tmp = (OP);             \
        Vec3_SET(DEST, _tmp, _tmp, _tmp);       \
} while (0)

#define Vec2_CLEAR(DEST) Vec2_SET(DEST,0,0)
#define Vec3_CLEAR(DEST) Vec3_SET(DEST,0,0,0)

#define Vec2_COPY(DEST, SRC)                                                   \
        Vec2_SET(DEST,                                                         \
                (SRC).x, (SRC).y)

#define Vec3_COPY(DEST, SRC)                                                   \
        Vec3_SET(DEST,                                                         \
                (SRC).x, (SRC).y, (SRC).z)    

#define Vec2_NEG(DEST, SRC)                                                    \
        Vec2_SET(DEST,                                                         \
                -(SRC).x, -(SRC).y)

#define Vec2_NEGTO(DEST) Vec2_NEG(DEST, DEST)

#define Vec3_NEG(DEST, SRC)                                                    \
        Vec3_SET(DEST,                                                         \
                -(SRC).x, -(SRC).y, -(SRC).z)    

#define Vec3_NEGTO(DEST) Vec3_NEG(DEST, DEST)


#define Vec2_ADD(DEST, V1, V2)                                                 \
        Vec2_SET(DEST,                                                         \
                (V1).x + (V2).x,                                             \
                (V1).y + (V2).y)

#define Vec2_ADDTO(DEST, V2)                                                   \
        (DEST).x += (V2).x,                                                  \
        (DEST).y += (V2).y

#define Vec3_ADD(DEST, V1, V2)                                                 \
        Vec3_SET(DEST,                                                         \
                (V1).x + (V2).x,                                             \
                (V1).y + (V2).y,                                             \
                (V1).z + (V2).z)                    

#define Vec3_ADDTO(DEST, V2)                                                   \
        (DEST).x += (V2).x,                                                  \
        (DEST).y += (V2).y,                                                  \
        (DEST).z += (V2).z


#define Vec2_SUB(DEST, V1, V2)                                                 \
        Vec2_SET(DEST,                                                         \
                (V1).x - (V2).x,                                             \
                (V1).y - (V2).y)

#define Vec2_SUBTO(DEST, V2)                                                   \
        (DEST).x -= (V2).x,                                                  \
        (DEST).y -= (V2).y

#define Vec3_SUB(DEST, V1, V2)                                                 \
        Vec3_SET(DEST,                                                         \
                (V1).x - (V2).x,                                             \
                (V1).y - (V2).y,                                             \
                (V1).z - (V2).z)              

#define Vec3_SUBTO(DEST, V2)                                                   \
        (DEST).x -= (V2).x,                                                  \
        (DEST).y -= (V2).y,                                                  \
        (DEST).z -= (V2).z


#define Vec2_MUL(DEST, V1, CONST)                                              \
        Vec2_SET(DEST,                                                         \
                (V1).x * (CONST),                                             \
                (V1).y * (CONST))

#define Vec2_MULTO(DEST, CONST)                                                \
        (DEST).x *= (CONST),                                                  \
        (DEST).y *= (CONST)

#define Vec3_MUL(DEST, V1, CONST)                                              \
        Vec3_SET(DEST,                                                         \
                (V1).x * (CONST),                                             \
                (V1).y * (CONST),                                             \
                (V1).z * (CONST))

#define Vec3_MULTO(DEST, CONST)                                                \
        (DEST).x *= (CONST),                                                  \
        (DEST).y *= (CONST),                                                  \
        (DEST).z *= (CONST)


#define Vec2_DIV(DEST, V1, CONST)                                              \
        Vec2_SET(DEST,                                                         \
                (V1).x / (CONST),                                             \
                (V1).y / (CONST))

#define Vec2_DIVTO(DEST, CONST)                                                \
        (DEST).x /= (CONST),                                                  \
        (DEST).y /= (CONST)

#define Vec3_DIV(DEST, V1, CONST)                                              \
        Vec3_SET(DEST,                                                         \
                (V1).x / (CONST),                                             \
                (V1).y / (CONST),                                             \
                (V1).z / (CONST))

#define Vec3_DIVTO(DEST, CONST)                                                \
        (DEST).x /= (CONST),                                                  \
        (DEST).y /= (CONST),                                                  \
        (DEST).z /= (CONST)


#define Vec2_DOT(VEC1, VEC2)                                                   \
       ((VEC1).x * (VEC2).x +                                                \
        (VEC1).y * (VEC2).y)

#define Vec3_DOT(VEC1, VEC2)                                                   \
       ((VEC1).x * (VEC2).x +                                                \
        (VEC1).y * (VEC2).y +                                                \
        (VEC1).z * (VEC2).z)

#define Vec2_SQR(VEC) Vec2_DOT(VEC, VEC)
#define Vec3_SQR(VEC) Vec3_DOT(VEC, VEC)

#define Vec2_ABSLEN(VEC) Vec_SQRT_FUNC(Vec2_SQR(VEC)) 
#define Vec3_ABSLEN(VEC) Vec_SQRT_FUNC(Vec3_SQR(VEC)) 

#define Vec2_COS(VEC1, VEC2) (                                                 \
        Vec2_DOT(VEC1,VEC2) / Vec_SQRT_FUNC(Vec2_SQR(VEC1) * Vec2_SQR(VEC2)))
#define Vec3_COS(VEC1, VEC2) (                                                 \
        Vec3_DOT(VEC1,VEC2) / Vec_SQRT_FUNC(Vec3_SQR(VEC1) * Vec3_SQR(VEC2)))

#define Vec2_ANG(VEC1, VEC2) (Vec_ACOS_FUNC(Vec2_COS(VEC1,VEC2)))
#define Vec3_ANG(VEC1, VEC2) (Vec_ACOS_FUNC(Vec3_COS(VEC1,VEC2)))

#define Vec2_CROSS1(V1, V2)                                                    \
        ((V1).y * (V2).x - (V1).y * (V2).x) 
#define Vec2_CROSS2(V1, V2)                                                    \
        ((V1).x * (V2).y - (V1).x * (V2).y) 

#define Vec2_CROSS(DEST, V1, V2)                                               \
        Vec_SET(DEST,                                                          \
                Vec2_CROSS1(V1,V2),                                            \
                Vec2_CROSS2(V1,V2))


#define Vec3_CROSS1(V1, V2)                                                    \
        ((V1).y * (V2).z - (V1).z * (V2).y) 
#define Vec3_CROSS2(V1, V2)                                                    \
        ((V1).z * (V2).x - (V1).x * (V2).z) 
#define Vec3_CROSS3(V1, V2)                                                    \
        ((V1).x * (V2).y - (V1).y * (V2).x)    

#define Vec3_CROSS(DEST, V1, V2)                                               \
        Vec_SET(DEST,                                                          \
                Vec3_CROSS1(V1,V2),                                            \
                Vec3_CROSS2(V1,V2),                                            \
                Vec3_CROSS3(V1,V2));

#define Vec2_EQ(V1, V2) \
	((V1).x == (V2).x && (V1).y == (V2).y)

#define Vec3_EQ(V1, V2) \
	((V1).x == (V2).x && (V1).y == (V2).y && (V1).z == (V2).z)

#ifdef __cplusplus
}
#endif 

/* specific vetor type used in this project */
#define vec_t Vec3D

#endif /*VECTOR_H_ */
