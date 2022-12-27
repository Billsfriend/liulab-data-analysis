library(CytoExploreR)

# setup compensation control
compenSet <- cyto_setup('220522-comp',
                        gatingTemplate = '220522-comp/enTempalate.csv')

# Transform fluorescent channels - default logicle transformations
compenSet <- cyto_transform(compenSet)

# Gate Cells
cyto_gate_draw(compenSet,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A", "SSC-A"))

# Gate single Cells
cyto_gate_draw(compenSet,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A", "FSC-H"))

# Compute spillover matrix
spill <- cyto_spillover_compute(compenSet,
                                parent = "Single Cells",
                                spillover = "220522-comp/Spillover-Matrix.csv")

# Visualise uncompensated data
cyto_plot_compensation(compenSet,
                       parent = "Single Cells")

# Visualise compensated data
cyto_plot_compensation(compenSet,
                       parent = "Single Cells",
                       spillover = "220522-comp/Spillover-Matrix.csv",
                       compensate = TRUE)

# Load and annotate samples
gatingSet <- cyto_setup("220522-stain",
                 gatingTemplate = "220522-stain/gatingTemplate.csv")

# Apply compensation
gatingSet <- cyto_compensate(gatingSet,
                      spillover = "220522-comp/Spillover-Matrix.csv")

# Transform fluorescent channels - default logicle transformations
gatingSet <- cyto_transform(gatingSet)

# Gate Cells
cyto_gate_draw(gatingSet,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))
# Gate Single Cells
cyto_gate_draw(gatingSet,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))

# Extract unstained control 
unstain <- cyto_extract(compenSet, "Single Cells")[[2]]

# gate live cells
cyto_gate_draw(gatingSet,
               parent = 'Single Cells',
               alias = 'Live Cells',
               channels = '7-AAD',
               overlay = unstain)

# gate CD45+ cells
cyto_gate_draw(gatingSet,
               parent = 'Live Cells',
               alias = 'CD45+ Cells',
               channels = 'CD45',
               overlay = unstain)

# gate 2D 7AAD- CD45+ cells
cyto_gate_draw(gatingSet,
               parent = 'Single Cells',
               alias = 'Live CD45+ Cells',
               channels = c('CD45','7-AAD'),
               type = "rectangle",
               overlay = unstain)

# Gating Tree
cyto_plot_gating_tree(gatingSet,
                      stat = "freq")

# Gating scheme
cyto_plot_gating_scheme(gatingSet,
                        back_gate = TRUE,
                        gate_track = TRUE)

# Compute medFI to a tibble
cyto_stats_compute(gatingSet[1:2],
                   alias = c("Live Cells",
                             "CD45+ Cells"),
                   stat = "median",
                   channels = c("CD45", "7-AAD"))
