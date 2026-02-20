// =============================================================================
// SCRIPT GMSH CORRIGIDO - MUFFLER TRANSFINITO (FLUXO CONTÍNUO)
// =============================================================================
Delete All; // Limpa geometria anterior para evitar sobreposição visual

// Parâmetros 0.005;
lc = 0.005;

H_tube   = 0.10; 
H_cavity = 0.20; 
H_total  = H_tube + 2*H_cavity; 

L_in      = 0.4;  
L_chamber = 0.6;  
L_gap     = 0.2;  
L_end     = 0.2;  
L_total   = L_in + L_chamber + L_gap + L_end;

// Coordenadas
Y0 = 0.0;
Y1 = H_cavity;
Y2 = H_cavity + H_tube;
Y3 = H_total;

X0 = 0.0;
X1 = L_in;
X2 = L_in + L_chamber;
X3 = X2 + L_gap;     
X4 = L_total;        

// Nós
n_in      = Ceil(L_in / lc) + 1;
n_chamber = Ceil(L_chamber / lc) + 1;
n_gap     = Ceil(L_gap / lc) + 1;
n_end     = Ceil(L_end / lc) + 1;

n_H_cav   = Ceil(H_cavity / lc) + 1;
n_H_tube  = Ceil(H_tube / lc) + 1;

// --- PONTOS ---
p1 = newp; Point(p1) = {X0, Y1, 0, lc}; 
p2 = newp; Point(p2) = {X1, Y1, 0, lc}; 
p3 = newp; Point(p3) = {X1, Y2, 0, lc}; 
p4 = newp; Point(p4) = {X0, Y2, 0, lc}; 

p5 = newp; Point(p5) = {X1, Y0, 0, lc}; 
p6 = newp; Point(p6) = {X2, Y0, 0, lc}; 
p7 = newp; Point(p7) = {X2, Y1, 0, lc}; 

p8 = newp; Point(p8) = {X2, Y2, 0, lc}; 
p9 = newp; Point(p9) = {X2, Y3, 0, lc}; 
p10= newp; Point(p10)= {X1, Y3, 0, lc}; 

p11= newp; Point(p11)= {X3, Y1, 0, lc}; 
p12= newp; Point(p12)= {X3, Y2, 0, lc}; 

p13= newp; Point(p13)= {X4, Y1, 0, lc}; 
p14= newp; Point(p14)= {X4, Y2, 0, lc}; 

// --- LINHAS ---
// Inlet
l_in_bot = newl; Line(l_in_bot) = {p1, p2};
l_in_R   = newl; Line(l_in_R)   = {p2, p3}; 
l_in_top = newl; Line(l_in_top) = {p3, p4};
l_in_L   = newl; Line(l_in_L)   = {p4, p1}; 

// Chamber Bot
l_cb_L   = newl; Line(l_cb_L)   = {p2, p5};
l_cb_bot = newl; Line(l_cb_bot) = {p5, p6};
l_cb_R   = newl; Line(l_cb_R)   = {p6, p7};
l_ch_floor=newl; Line(l_ch_floor)={p2, p7}; 

// Channel Vert Interface
l_ch_out = newl; Line(l_ch_out) = {p7, p8}; 

// Chamber Top
l_ct_floor=newl; Line(l_ct_floor)={p8, p3}; // Definida da Dir para Esq
l_ct_L   = newl; Line(l_ct_L)   = {p3, p10};
l_ct_top = newl; Line(l_ct_top) = {p10, p9};
l_ct_R   = newl; Line(l_ct_R)   = {p9, p8};

// Outlets
l_o1_bot = newl; Line(l_o1_bot) = {p7, p11};
l_target = newl; Line(l_target) = {p11, p12}; 
l_o1_top = newl; Line(l_o1_top) = {p12, p8};

l_o2_bot = newl; Line(l_o2_bot) = {p11, p13};
l_end    = newl; Line(l_end)    = {p13, p14}; 
l_o2_top = newl; Line(l_o2_top) = {p14, p12};


// --- LOOPS E SUPERFÍCIES ---

// Inlet (CCW)
cl_in = newll; Curve Loop(cl_in) = {l_in_bot, l_in_R, l_in_top, l_in_L};
s_in  = news;  Plane Surface(s_in) = {cl_in};

// Cavity Bot (CCW)
cl_cb = newll; Curve Loop(cl_cb) = {l_cb_bot, l_cb_R, -l_ch_floor, l_cb_L};
s_cb  = news;  Plane Surface(s_cb) = {cl_cb};

// Channel (CCW)
cl_ch = newll; Curve Loop(cl_ch) = {l_ch_floor, l_ch_out, l_ct_floor, -l_in_R};
s_ch  = news;  Plane Surface(s_ch) = {cl_ch};

// Cavity Top (CORRIGIDO PARA CCW)
// Caminho: p3 -> p8 -> p9 -> p10 -> p3
// l_ct_floor é p8->p3. Precisamos p3->p8 (Inverso: -l_ct_floor)
// l_ct_R é p9->p8. Precisamos p8->p9 (Inverso: -l_ct_R)
// l_ct_top é p10->p9. Precisamos p9->p10 (Inverso: -l_ct_top)
// l_ct_L é p3->p10. Precisamos p10->p3 (Inverso: -l_ct_L)
cl_ct = newll; Curve Loop(cl_ct) = {-l_ct_floor, -l_ct_R, -l_ct_top, -l_ct_L};
s_ct  = news;  Plane Surface(s_ct) = {cl_ct};

// Outlet 1 (CCW)
cl_o1 = newll; Curve Loop(cl_o1) = {l_o1_bot, l_target, l_o1_top, -l_ch_out};
s_o1  = news;  Plane Surface(s_o1) = {cl_o1};

// Outlet 2 (CCW)
cl_o2 = newll; Curve Loop(cl_o2) = {l_o2_bot, l_end, l_o2_top, -l_target};
s_o2  = news;  Plane Surface(s_o2) = {cl_o2};


// =============================================================================
// 4. MALHA TRANSFINITA
// =============================================================================

Transfinite Curve {l_in_bot, l_in_top} = n_in;
Transfinite Curve {l_cb_bot, l_ch_floor, l_ct_floor, l_ct_top} = n_chamber;
Transfinite Curve {l_o1_bot, l_o1_top} = n_gap;
Transfinite Curve {l_o2_bot, l_o2_top} = n_end;

Transfinite Curve {l_in_L, l_in_R, l_ch_out, l_target, l_end} = n_H_tube;
Transfinite Curve {l_cb_L, l_cb_R, l_ct_L, l_ct_R} = n_H_cav;

// Superfícies - CORRIGIDAS PARA O NOVO LOOP
// Ordem dos cantos deve seguir o loop definido acima

Transfinite Surface {s_in} = {p1, p2, p3, p4};
Transfinite Surface {s_cb} = {p5, p6, p7, p2};
Transfinite Surface {s_ch} = {p2, p7, p8, p3};

// CORREÇÃO: Começa em p3 e segue CCW (p3->p8->p9->p10)
Transfinite Surface {s_ct} = {p3, p8, p9, p10}; 

Transfinite Surface {s_o1} = {p7, p11, p12, p8};
Transfinite Surface {s_o2} = {p11, p13, p14, p12};

Recombine Surface {s_in, s_cb, s_ch, s_ct, s_o1, s_o2};

// =============================================================================
// 5. GRUPOS FÍSICOS
// =============================================================================

Physical Surface("Material,Ar,1,1.21,343.0,0.0") = {s_in, s_cb, s_ch, s_ct, s_o1, s_o2};
Physical Surface("Fixed,0.0") = {s_in, s_ch, s_o1, s_o2};

Physical Curve("Vn,1E-3,0.0,0.0") = {l_in_L};       
Physical Curve("Pressure,0.0,0.0,0.0") = {l_end}; 
Physical Curve("Target") = {l_target};               

Mesh.Algorithm = 8; 
Mesh 2;
Save "muffler.msh";