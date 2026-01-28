SetFactory("OpenCASCADE");

// Tamanho do elemento (lc)
// Para 1000 Hz, o comprimento de onda é ~0.34m. 
// lc = 0.005m garante > 60 elementos por comprimento de onda (muito preciso).
lc = 0.005;

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// Dimensões do Duto
H = 0.2;       // Altura do duto (m)
L_in = 0.3;    // Fim da região de entrada (m)
L_d  = 0.7;    // Fim da região de projeto (m)
L_tgt = 0.9;   // Posição da linha de medição (Target) (m)
L_out = 1.0;   // Fim do duto (Outlet) (m)

// -------------------------------------------------------------------
// Pontos
// -------------------------------------------------------------------
// Base (y=0)
Point(1) = {0,     0, 0, lc};
Point(2) = {L_in,  0, 0, lc};
Point(3) = {L_d,   0, 0, lc};
Point(4) = {L_tgt, 0, 0, lc};
Point(5) = {L_out, 0, 0, lc};

// Topo (y=H)
Point(6)  = {0,     H, 0, lc};
Point(7)  = {L_in,  H, 0, lc};
Point(8)  = {L_d,   H, 0, lc};
Point(9)  = {L_tgt, H, 0, lc};
Point(10) = {L_out, H, 0, lc};

// -------------------------------------------------------------------
// Linhas
// -------------------------------------------------------------------
// Horizontais (Base)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};

// Horizontais (Topo)
Line(5) = {6, 7};
Line(6) = {7, 8};
Line(7) = {8, 9};
Line(8) = {9, 10};

// Verticais
Line(9)  = {1, 6};   // Entrada (Inlet)
Line(10) = {2, 7};   // Início do Domínio de Projeto
Line(11) = {3, 8};   // Fim do Domínio de Projeto
Line(12) = {4, 9};   // *** LINHA ALVO (TARGET) ***
Line(13) = {5, 10};  // Saída (Outlet)

// -------------------------------------------------------------------
// Superfícies
// -------------------------------------------------------------------
Curve Loop(100) = {1, 10, -5, -9}; // Região de Entrada (Ar)
Plane Surface(1) = {100};

Curve Loop(200) = {2, 11, -6, -10}; // Região de PROJETO
Plane Surface(2) = {200};

Curve Loop(300) = {3, 12, -7, -11}; // Região de Saída 1 (Ar)
Plane Surface(3) = {300};

Curve Loop(400) = {4, 13, -8, -12}; // Região de Saída 2 (Ar)
Plane Surface(4) = {400};

// Garante que os nós nas interfaces (ex: entre surf 1 e 2) sejam os mesmos
Coherence;

// -------------------------------------------------------------------
// Definições Físicas para o Código Julia
// -------------------------------------------------------------------

// 1. Materiais (Tag, Nome, ID, Densidade, c, mu)
// O Ar preenche todo o domínio
Physical Surface("Material,Ar,1,1.21,343.0,0.0") = {1, 2, 3, 4};

// O Sólido só precisa existir como opção na região de projeto (Surf 2)
Physical Surface("Material,Solido,2,2700.0,5000.0,0.0") = {2};

// 2. Elementos Fixos (Não-projeto)
// Superfícies 1, 3 e 4 são travadas como Ar (gama = 0)
Physical Surface("Fixed,0.0") = {1, 3, 4};

// 3. Condições de Contorno
// Excitação de Velocidade Normal na esquerda (Amplitude 1.0 m/s)
Physical Curve("Vn,1.0,0.0,0.0") = {9};

// Pressão Nula (Release) na direita
Physical Curve("Pressure,0.0,0.0,0.0") = {13};

// As paredes superior e inferior ficam rígidas automaticamente (Neumann homogênea do MEF)

// 4. Nós Alvo para calcular o SPL
// Mede a pressão média que cruza a linha x=0.9m
Physical Curve("Target") = {12};


// -------------------------------------------------------------------
// Malha (Quad)
// -------------------------------------------------------------------
Recombine Surface {1, 2, 3, 4};
Mesh.Algorithm = 8; // Frontal-Delaunay for quads

Mesh 2;
Save "duct_example.msh";