SetFactory("OpenCASCADE");

// Tamanho do elemento (lc)
// 0.005m garante excelente resolução para os efeitos de difração no canto
lc = 0.005;

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// -------------------------------------------------------------------
// Parâmetros Geométricos (em metros)
// -------------------------------------------------------------------
W = 0.2; // Largura do duto

// Coordenadas X
X0 = 0.0;       // Início da entrada
X1 = 0.4;       // Início do canto
X2 = X1 + W;    // Fim do canto (0.6)
X3 = X2 + W;    // Fim da cavidade exterior (0.8)

// Coordenadas Y
Y0 = 0.0;       // Fundo da cavidade exterior
Y1 = W;         // Início do canto (0.2)
Y2 = Y1 + W;    // Fim do canto (0.4)
Y3 = 0.7;       // Linha Alvo (Target)
Y4 = 0.8;       // Saída (Outlet)

// -------------------------------------------------------------------
// Pontos
// -------------------------------------------------------------------
// Nível Y0 (Fundo)
Point(1) = {X1, Y0, 0, lc};
Point(2) = {X2, Y0, 0, lc};
Point(3) = {X3, Y0, 0, lc};

// Nível Y1 
Point(4) = {X1, Y1, 0, lc};
Point(5) = {X2, Y1, 0, lc};
Point(6) = {X3, Y1, 0, lc};

// Nível Y2 (Entrada)
Point(7) = {X0, Y2, 0, lc};
Point(8) = {X1, Y2, 0, lc};
Point(9) = {X2, Y2, 0, lc};
Point(10)= {X3, Y2, 0, lc};

// Nível Y2_bottom (Entrada)
Point(11) = {X0, Y1, 0, lc};

// Nível Y3 (Target)
Point(12) = {X1, Y3, 0, lc};
Point(13) = {X2, Y3, 0, lc};

// Nível Y4 (Outlet)
Point(14) = {X1, Y4, 0, lc};
Point(15) = {X2, Y4, 0, lc};


// -------------------------------------------------------------------
// Linhas
// -------------------------------------------------------------------
// Horizontais
Line(1) = {1, 2}; // Fundo Cavidade 1
Line(2) = {2, 3}; // Fundo Cavidade 2

Line(3) = {11, 4}; // Fundo Entrada
Line(4) = {4, 5};  // Divisa Cavidade/Canto
Line(5) = {5, 6};  // Fundo Cavidade 3

Line(6) = {7, 8};  // Topo Entrada
Line(7) = {8, 9};  // Divisa Canto/Saída
Line(8) = {9, 10}; // Topo Cavidade 3

Line(9)  = {12, 13}; // *** TARGET LINE ***
Line(10) = {14, 15}; // Outlet Line

// Verticais
Line(11) = {11, 7}; // Inlet
Line(12) = {1, 4};  // Esquerda Cavidade 1
Line(13) = {4, 8};  // Esquerda Canto
Line(14) = {8, 12}; // Esquerda Saída 1
Line(15) = {12, 14}; // Esquerda Saída 2

Line(16) = {2, 5};  // Direita Cavidade 1 / Esq Cavidade 2
Line(17) = {5, 9};  // Direita Canto / Esq Cavidade 3
Line(18) = {9, 13}; // Direita Saída 1
Line(19) = {13, 15}; // Direita Saída 2

Line(20) = {3, 6};  // Direita Cavidade 2
Line(21) = {6, 10}; // Direita Cavidade 3


// -------------------------------------------------------------------
// Superfícies (Curve Loops)
// -------------------------------------------------------------------
// Duto Principal (Ar Fixo)
Curve Loop(100) = {3, 13, -6, -11}; Plane Surface(1) = {100}; // Duto de Entrada
Curve Loop(200) = {4, 17, -7, -13}; Plane Surface(2) = {200}; // Canto do Fluido
Curve Loop(300) = {7, 18, -9, -14}; Plane Surface(3) = {300}; // Duto Saída 1
Curve Loop(400) = {9, 19, -10, -15}; Plane Surface(4) = {400}; // Duto Saída 2

// Domínio de Projeto (Cavidade Exterior em "L")
Curve Loop(500) = {1, 16, -4, -12}; Plane Surface(5) = {500}; // Cavidade Baixo
Curve Loop(600) = {2, 20, -5, -16}; Plane Surface(6) = {600}; // Cavidade Canto (Extra)
Curve Loop(700) = {5, 21, -8, -17}; Plane Surface(7) = {700}; // Cavidade Direita

Coherence;

// -------------------------------------------------------------------
// Grupos Físicos para o Código Julia
// -------------------------------------------------------------------

// 1. Materiais
// Todo o domínio começa como Ar
Physical Surface("Material,Ar,1,1.21,343.0,0.0") = {1, 2, 3, 4, 5, 6, 7};

// 2. Elementos Fixos (Não-projeto)
// Superfícies 1, 2, 3 e 4 são o caminho do duto: Travadas como Ar
Physical Surface("Fixed,0.0") = {1, 2, 3, 4};

// 3. Condições de Contorno
// Excitação na entrada (v_n = 1.0 m/s)
Physical Curve("Vn,1.0,0.0,0.0") = {11};

// Pressão Nula na saída superior
Physical Curve("Pressure,0.0,0.0,0.0") = {10};

// 4. Nós Alvo
// Mede a pressão que passa pela linha y = 0.7
Physical Curve("Target") = {9};

// -------------------------------------------------------------------
// Malha (Quad)
// -------------------------------------------------------------------
Recombine Surface {1, 2, 3, 4, 5, 6, 7};
Mesh.Algorithm = 8; 

Mesh 2;
Save "bend_cavity.msh";