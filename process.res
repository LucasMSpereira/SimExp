InOut.Extract

mtx.formnf

nf:
Int32[0, 1, 2, 0, 4, 6, 0, 9, 11]
Int32[0, 0, 0, 3, 5, 7, 8, 10, 12]
Number of equations = 12

coords:
Float32[-0.5, -0.5]
Float32[0.0, -0.5]
Float32[0.5, -0.5]
Float32[-0.5, 0.0]
Float32[0.0, 0.0]
Float32[0.5, 0.0]
Float32[-0.5, 0.5]
Float32[0.0, 0.5]
Float32[0.5, 0.5]

g_num:
Int32[1, 1, 2, 2, 4, 4, 5, 5]
Int32[2, 5, 3, 6, 5, 8, 6, 9]
Int32[5, 4, 6, 5, 8, 7, 9, 8]

g_g:
Int32[0, 0, 1, 1, 0, 0, 4, 4]
Int32[0, 0, 0, 0, 3, 3, 5, 5]
Int32[1, 4, 2, 6, 4, 9, 6, 11]
Int32[0, 5, 0, 7, 5, 10, 7, 12]
Int32[4, 0, 6, 4, 9, 0, 11, 9]
Int32[5, 3, 7, 5, 10, 8, 12, 10]

/////////////////////////////////////////////

     ELEMENTS STIFFNESS

/////////////////////////////////////////////

num.sample

points:
[0.3333333333333333, 0.3333333333333333]
weights:
[0.5]

num.stiffSerialAssembly

mtx.fkdiag

kdiag = Int32[1, 3, 4, 8, 13, 19, 26, 32, 39, 47, 55, 64]

num.deemat

dee:
Float32[1.346154f6, 576923.2, 0.0]
Float32[576923.2, 1.346154f6, 0.0]
Float32[0.0, 0.0, 384615.4]

/////////////////  SITFFNESS OF ELEMENT 1 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 2 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 3 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 4 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 5 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 6 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 7 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

num.fsparv


/////////////////  SITFFNESS OF ELEMENT 8 //////////////

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

num.fsparv


num.sparin


num.spabac

Displacements:
Float32[1.9499998f-7]
Float32[3.8999997f-7]
Float32[-4.5499993f-7]
Float32[1.95f-7]
Float32[-4.5499993f-7]
Float32[3.9f-7]
Float32[-4.5499993f-7]
Float32[-9.0999987f-7]
Float32[1.9500001f-7]
Float32[-9.0999987f-7]
Float32[3.9000003f-7]
Float32[-9.1f-7]

InOut.DispOut

/////////////////////////////////////////////

Calculating stresses at integrating points

/////////////////////////////////////////////

num.deemat

dee:
Float32[1.346154f6, 576923.2, 0.0]
Float32[576923.2, 1.346154f6, 0.0]
Float32[0.0, 0.0, 384615.4]

InOut.HeaderStress

Integrating point 1 of element 1

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 2

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 3

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 4

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 5

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 6

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 7

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[0.0, 0.5]
Float32[-0.5, 0.0]

bee:
[-2.0, 0.0, 2.0, 0.0, 0.0, 0.0]
[0.0, 0.0, 0.0, -2.0, 0.0, 2.0]
[0.0, -2.0, -2.0, 2.0, 2.0, 0.0]

Jacobian determinant = 0.25

InOut.StressOut

Integrating point 1 of element 8

num.shapeDer

der:
[0, -1, 1]
[1, -1, 0]

num.beemat

Jacobian:
Float32[-0.5, 0.0]
Float32[-0.5, -0.5]

bee:
[0.0, 0.0, 2.0, 0.0, -2.0, 0.0]
[0.0, -2.0, 0.0, 0.0, 0.0, 2.0]
[-2.0, 0.0, 0.0, 2.0, 2.0, -2.0]

Jacobian determinant = 0.25

InOut.StressOut
