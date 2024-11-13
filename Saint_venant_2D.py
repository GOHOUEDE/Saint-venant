import numpy as np
2 import math
3 import matplotlib . pyplot as plt
4 t = np . array ([0 ,0.2 ,0.4 ,0.6 ,0.8 ,1 ,1.2 ,1.4 ,1.6 ,1.8])
5 il = np . arange (0 ,10 ,1)
6
7 D =5
8 io =0.1
9 ho =1.5
10 nu =0.03
11 g =9.8
12 Ho = ho
13
14 for i in range (len( t ) ) :
15 beta =1 -(2* ho / D )
16 bo = D * np . sqrt (1 - beta **2)
17 R = D /4*(1 -( beta * np . sqrt (1 - beta **2) / math . acos ( beta ) ) )
18 Cw = R **(1/6) / nu
19 x =10+(12* ho / bo ) /(3+(6* ho / bo ) )
20 vo = Cw * np . sqrt ( io * R )
21 z = Cw * t + il
22 a = Cw - vo
23 b =2* Cw -( x * vo )
24 A =2*( a * Ho ) -( ho * vo ) / b
25 B = Ho - A
26 c = -1* g * io * b /( vo * g * ho - a **2)
27 C = np . exp ( c * z )
28 H = A +( B * C )
29
30 V =1/ ho * a *( H - Ho ) + vo
31 plt . plot ( il ,H , label = f’courbe de H pour t={t[i]} ’)
32 plt . plot ( il ,V , label = f’courbe de V pour t=={t[i]} ’)
33 plt . xlabel (’Axe de x’)
34 plt . ylabel (’ H et V’)
35 plt . legend ()
36 plt . show ()
