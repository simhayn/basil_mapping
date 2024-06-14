####linkage map view####
# Step 1: Install and load the LinkageMapView package if not already installed
if (!requireNamespace("LinkageMapView", quietly = TRUE)) {
  install.packages("LinkageMapView")
}
library(LinkageMapView)

# Load other necessary packages
library(qtl)


# Assuming 'qtldf.aq' is your QTL summary dataframe as prepared in previous steps
# Identify the linkage groups that have QTLs for the traits of interest
traits_of_interest <- names(qtlist)[1:6]  # Replace with your traits
# Find the LGs with QTLs
(lg_with_qtls <- unique(unlist(lapply(qtlist, function(x) if(!is.null(x)) x$chr))))
(lg_with_qtls<-lg_with_qtls[1:nrow(qtlist$AFF)])

# Step 3: Filter the cross object to include only linkage groups with QTLs
filtered_cross <- subset(data, chr = lg_with_qtls)

# Step 4: Plot and save to PDF

output_file <- file.path("results",file= "linkage_map_with_qtls.pdf")

# Generate the linkage map plot with QTLs
lmv.linkage.plot(data,output_file,mapthese=lg_with_qtls,qtlscanone = scan1,dupnbr = T)

# Close the PDF device
dev.off()

