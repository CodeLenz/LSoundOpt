// 1. Ocultar as vistas originais
View[0].Visible = 0;
View[1].Visible = 0;

// 2. Executar o filtro na vista elementar (View 0)
Plugin(ExtractElements).View = 0;
Plugin(ExtractElements).MinVal = 1;
Plugin(ExtractElements).MaxVal = 1;
Plugin(ExtractElements).Run;

// Capturar o índice da vista extraida (elementar)
view_extraida = PostProcessing.NbViews - 1;
View[view_extraida].Visible = 0; // Oculta ela também

// 3. Converter a vista extraida de Elementar para Nodal
Plugin(Smooth).View = view_extraida;
Plugin(Smooth).Run;

// Capturar o índice da nova vista suavizada (nodal)
view_suavizada = PostProcessing.NbViews - 1;

// 4. Forçar a nova vista nodal a buscar as cores no campo original (View 1)
View[view_suavizada].ExternalView = 1;