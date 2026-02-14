SetFactory("OpenCASCADE");

// Element size
// 0.005m ensures >60 elements per wavelength at 1000 Hz.
lc = 0.005;

Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// Geometry Parameters
H_domain   = 0.10; // Height of Top/Bottom Design Domains
H_channel  = 0.10; // Height of Central Air Channel
H_total    = 2*H_domain + H_channel; // Total height (0.3 m)

L_in  = 0.3;      // Inlet region length
L_d   = 0.7;      // Design domain end
L_tgt = 0.9;      // Target measurement line
L_out = 1.0;      // Outlet

// Y-coordinates for layers
Y0 = 0.0;
Y1 = H_domain;              // Bottom of central channel
Y2 = H_domain + H_channel;  // Top of central channel
Y3 = H_total;               // Top of duct

// -------------------------------------------------------------------
// Points
// -------------------------------------------------------------------
// Bottom (y=Y0)
Point(1) = {0, Y0, 0, lc}; Point(2) = {L_in, Y0, 0, lc}; Point(3) = {L_d, Y0, 0, lc}; Point(4) = {L_tgt, Y0, 0, lc}; Point(5) = {L_out, Y0, 0, lc};

// Channel Bottom (y=Y1)
Point(6) = {0, Y1, 0, lc}; Point(7) = {L_in, Y1, 0, lc}; Point(8) = {L_d, Y1, 0, lc}; Point(9) = {L_tgt, Y1, 0, lc}; Point(10) = {L_out, Y1, 0, lc};

// Channel Top (y=Y2)
Point(11) = {0, Y2, 0, lc}; Point(12) = {L_in, Y2, 0, lc}; Point(13) = {L_d, Y2, 0, lc}; Point(14) = {L_tgt, Y2, 0, lc}; Point(15) = {L_out, Y2, 0, lc};

// Top (y=Y3)
Point(16) = {0, Y3, 0, lc}; Point(17) = {L_in, Y3, 0, lc}; Point(18) = {L_d, Y3, 0, lc}; Point(19) = {L_tgt, Y3, 0, lc}; Point(20) = {L_out, Y3, 0, lc};

// -------------------------------------------------------------------
// Lines
// -------------------------------------------------------------------
// Horizontal Lines
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 5};       // Y0
Line(5) = {6, 7}; Line(6) = {7, 8}; Line(7) = {8, 9}; Line(8) = {9, 10};      // Y1
Line(9) = {11, 12}; Line(10) = {12, 13}; Line(11) = {13, 14}; Line(12) = {14, 15}; // Y2
Line(13) = {16, 17}; Line(14) = {17, 18}; Line(15) = {18, 19}; Line(16) = {19, 20}; // Y3

// Vertical Lines
Line(17) = {1, 6}; Line(18) = {6, 11}; Line(19) = {11, 16};   // Inlet
Line(20) = {2, 7}; Line(21) = {7, 12}; Line(22) = {12, 17};   // Start Design
Line(23) = {3, 8}; Line(24) = {8, 13}; Line(25) = {13, 18};   // End Design
Line(26) = {4, 9}; Line(27) = {9, 14}; Line(28) = {14, 19};   // Target
Line(29) = {5, 10}; Line(30) = {10, 15}; Line(31) = {15, 20}; // Outlet

// -------------------------------------------------------------------
// Surfaces (Curve Loops)
// -------------------------------------------------------------------
// --- INLET REGION ---
Curve Loop(101) = {1, 20, -5, -17}; Plane Surface(1) = {101}; // Bottom
Curve Loop(102) = {5, 21, -9, -18}; Plane Surface(2) = {102}; // Channel
Curve Loop(103) = {9, 22, -13, -19}; Plane Surface(3) = {103}; // Top

// --- CENTRAL REGION (DESIGN + CHANNEL) ---
Curve Loop(201) = {2, 23, -6, -20}; Plane Surface(4) = {201}; // BOTTOM DESIGN DOMAIN
Curve Loop(202) = {6, 24, -10, -21}; Plane Surface(5) = {202}; // CENTRAL CHANNEL (FIXED)
Curve Loop(203) = {10, 25, -14, -22}; Plane Surface(6) = {203}; // TOP DESIGN DOMAIN

// --- OUTLET REGION 1 ---
Curve Loop(301) = {3, 26, -7, -23}; Plane Surface(7) = {301}; // Bottom
Curve Loop(302) = {7, 27, -11, -24}; Plane Surface(8) = {302}; // Channel
Curve Loop(303) = {11, 28, -15, -25}; Plane Surface(9) = {303}; // Top

// --- OUTLET REGION 2 ---
Curve Loop(401) = {4, 29, -8, -26}; Plane Surface(10) = {401}; // Bottom
Curve Loop(402) = {8, 30, -12, -27}; Plane Surface(11) = {402}; // Channel
Curve Loop(403) = {12, 31, -16, -28}; Plane Surface(12) = {403}; // Top

Coherence;

// -------------------------------------------------------------------
// Physical Groups
// -------------------------------------------------------------------

// 1. All Air by default
Physical Surface("Material,Ar,1,1.21,343.0,0.0") = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

// 3. Fixed Elements (Everything EXCEPT 4 and 6)
Physical Surface("Fixed,0.0") = {1, 2, 3, 5, 7, 8, 9, 10, 11, 12};

// 4. Boundary Conditions
// Inlet Source
Physical Curve("Vn,1E-3,0.0,0.0") = {17, 18, 19};

// Outlet (Pressure Release)
Physical Curve("Pressure,0.0,0.0,0.0") = {29, 30, 31};

// 5. Target Nodes for SPL (x=0.9, across full height)
Physical Curve("Target") = {26, 27, 28};

// -------------------------------------------------------------------
// Mesh
// -------------------------------------------------------------------
Recombine Surface {1:12};
Mesh.Algorithm = 8; 

Mesh 2;
Save "duct_mirrored.msh";