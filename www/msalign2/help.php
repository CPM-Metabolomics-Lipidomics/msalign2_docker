<?PHP
include("header.inc.php");
echo $beginheader;
echo $scriptheader;
echo $bodyheader;

echo "<p><b><u>How to use :</u></b>\n";
echo "<ul>\n";
echo "<li>Before you can start your alignment your files need to be in the correct location. The website can only locate your files if they are located in <b>\\\paradata\work\webtools\msalign</b>. This location is linked to the folder <b>data</b> seen on the homepage.</li>\n";
echo "<li>First select a <b>Master file</b> to which you want to align all other files. It is only possible to select one file as a Master file.</li>\n";
echo "<li>Next, select all files you want to use fo alignment.</li>\n";
echo "<li>Select the parameters you want to use. The default parameters already set should give already a proper alignment.\n";
echo "<ul>\n";
echo "<li><b>Mass measurement error [ppm]:</b> Set the mass error of you measurements in ppm</li>\n";
echo "<li><b>Standard deviation retention time [scans]:</b> Set the standard deviation of the retention time of your measure in number of scans.</li>\n";
echo "<li><b>Number of features :</b> Set the number of features you want to use for aligning. You need at least 100 features for alignment. The more complex your samples is the more features you want to use.</li>\n";
echo "<li><b>Start at scan number :</b> From which scan number to start the alignment.</li>\n";
echo "<li><b>End at scan number :</b> Until which scan number to end the alignment.</li>\n";
echo "<li><b>Endpoint X :</b> The alignment starts at (0,0), but also needs an endpoint beyond the last scan, this is the X value. e.g. add 500 - 1000 to the end scan.</li>\n";
echo "<li><b>Endpoint Y :</b> The alignment starts at (0,0), but also needs an endpoint beyond the last scan, this is the Y value. e.g. add 500 - 1000 to the end scan.</li>\n";
echo "<li><b>Costs per breakpoint [0 - 1]:</b> Set the costs per breakpoint. This determines if a breakpoint is kept or discarded. Lower this value to keep more breaktpoint. Default is 0.5 and should generally give a good aligment.</li>\n";
echo "<li><b>Fitness [1 - ..]:</b> The fitness determines ... Increasing this value also increases the possibility that an outlier is used for alignment. Generally the default value 1 should give good alignment.</li>\n";
echo "</ul>\n";
echo "</li>\n";
echo "<li>After starting the alignment, you can refresh the statuspage to monitor the progress. This page can also be bookmarked!</li>\n";
echo "<li>After a file is aligned a link appears which will give you a link to show you how good / bad your alignment is. <a href=\"./images/good_example.png\">Good</a> and a <a href=\"./images/bad_example.png\">bad</a> aligment example</li>\n";
echo "<li>After the alignment is finished a button for removing the temporary files appears. Please use this button to remove the temporary files and be redirected to the homepage again.</li>\n";
echo "<li>Your aligned mzXML files will be located in the same folder as your original files.</li>\n";
echo "</ul>\n";
echo "</p>\n";

include("footer.inc.php");
echo $footer;
?>
