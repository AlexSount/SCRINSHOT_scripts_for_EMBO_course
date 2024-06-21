// set the input (correspond to the results of "5_count_matrix_creation_v3.R"script) and output directories 
// output direcrory should not be in the input
// in the line-9, set the x,y dimensions of the analyzed dataset (in pixels) 
// in the line-25, set the path of the ROI Manager *.zip file , with the cell-ROIs (not the nuclear)
input = "C:/Users/alex/Desktop/SCRINSHOT_updated_tutorial/analysis_20220503/output/single_gene_dot_coordinates/";
output = "C:/Users/alex/Desktop/SCRINSHOT_updated_tutorial/analysis_20220503/output/roi_counts/";

function action(input, output, filename) {
        newImage("blank", "8-bit black", 6700, 5000, 1);
        open(input + filename);
        Table.rename(filename, "Results");
		numberOfRows = nResults;
        for (i=0; i<numberOfRows; ++i) {
		x = getResult("Location_Center_X", i);
		y = getResult("Location_Center_Y", i);
		run("Specify...", "width=1 height=1 x=x y=y");
		setBackgroundColor(255, 255, 255);
		run("Clear", "slice");
		setForegroundColor(255, 255, 255);
		run("Draw", "slice");
		}
	selectWindow("Results"); 
    run("Close" );
    
	roiManager("Open", "C:/Users/alex/Desktop/SCRINSHOT/SCRINSHOT_paper/analysis_procedure/507_s7/507_s7_all_cell_rois.zip");
	roiManager("Measure");
	cells = roiManager("count");
	array = newArray(cells);
 	 for (i=0; i<array.length; i++) {
      array[i] = i;
	  }
		
	roiManager("Select", array);
	roiManager("Delete");
	saveAs("Results", output + filename);    
    selectWindow("Results"); 
    run("Close" );
        close();
}

setBatchMode(true); 
list = getFileList(input);
for (j = 0; j < list.length; j++)
        action(input, output, list[j]);
setBatchMode(false);
