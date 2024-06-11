<?PHP
//include("functions.inc.php");    // include several functions


//---------------------------------------------------------------------------
if ((!isset($_GET['Settings'])) && (!isset($_POST['Settings'])))
{
    $Settings[Status] = "start";
}
else
{
  $Settings = unserialize(base64_decode($_GET['Settings']));
}
if (isset($_POST['Settings']))
{
  $Settings = unserialize(base64_decode($_POST['Settings']));
  if (isset($_POST['MasterFile']))
  {
   $Settings[MasterFile] = $_POST['MasterFile'];
  }
  if (isset($_POST['Files']))
  {
    $Settings[Files] = $_POST['Files'];
  }
}


  include("header.inc.php");
  echo $beginheader;
  echo $scriptheader;
  echo $bodyheader;

if ($Settings[Status] != "parameters")
{
  $DataDir = urldecode($Settings[dir]);
  if ($DataDir == "") { $DataDir = "./data"; }
  $ShowDirs = explode("/", $DataDir);

  echo "You are at : <b>&nbsp;/&nbsp;";
  unset($ShowDirs[0]); // remove the dot
  $folder = "./"; // initialize the folder
  foreach ($ShowDirs as $Directory)
  {
    $folder .= $Directory;
    echo "<span style=\"background-color: rgb(255,255,255);\" onmouseover=\"this.style.backgroundColor='#C4C4C4'\" onmouseout=\"this.style.backgroundColor='#FFFFFF'\" bgcolor=\"#FFFFFF\">";
    $Settings[dir] = urlencode($folder);
    echo "<a href=\"index.php?Settings=" .base64_encode(serialize($Settings)) ."\">" .$Directory ."</a></span>&nbsp;/&nbsp;";
    $folder .= "/";
  }
  echo "</b><br><br>";

  $ListDirs = "find ./data/ -maxdepth 1 -type d";
  exec ($ListDirs, $Dirs);
  sort($Dirs);
  $ListSubDirs = "find " .str_replace(" ","\ ",$DataDir) ." -maxdepth 1 -type d";
  exec ($ListSubDirs, $SubDirs);
  sort($SubDirs);
  
  if ($Settings[Status] == "start")
  {
    echo "<center><h3>Select a master file :</h3></center>";
    echo "<form method=\"post\" action=\"index.php\">";
  }
  if ($Settings[Status] == "master")
  {
    echo "<p><b>Selected masterfile :</b> " .$Settings[MasterFile] ."</p>";
    echo "<center><h3>Select files for aligning :</h3></center>";
    echo "<form method=\"post\" action=\"index.php\">";
  }
  echo "<hr>";
  echo "<table width=\"100%\">";
    echo "<tr valign=\"top\">";
    echo "<td width=\"50%\">";
      echo "<table>";
        if(count($SubDirs)>1)
	{
          unset($SubDirs[0]); // remove the first value
	}
        foreach ($SubDirs as $Value)
        {
	  if($Value==$DataDir)
	  {
	    echo "<tr><td>..</td></tr>";
	  }
	  else
	  {
	    echo "<tr style=\"background-color: rgb(255,255,255);\" onmouseover=\"this.style.backgroundColor='#C4C4C4'\" onmouseout=\"this.style.backgroundColor='#FFFFFF'\" bgcolor=\"#FFFFFF\"><td>";
	    $Settings[dir] = urlencode($Value);
            echo "<a href=\"index.php?Settings=" .base64_encode(serialize($Settings)) ."\">" .$Value ."</a>";
	    echo "</td></tr>";
	  }
        }
        echo "</table>";
      echo "</td><td width=\"50%\">";
        // formulate the command, make sure to add > for redirecting output to a file
        $ListFiles = "find " .str_replace(" ","\ ",$DataDir) ."/ -iname \"*.mzxml\"" ." -maxdepth 1 -type f";
        // executes the command and writes the input to a file
        exec($ListFiles, $FileNames);
        sort($FileNames);

        echo "<table>";
        echo "<tr><td><b>Please select file(s) to submit :</b>&nbsp;&nbsp;&nbsp;";
        echo "<input type=\"Submit\" name=\"Submit\" value=\"Submit\"></td></tr>";

        $cbCounter = 0;
        foreach ($FileNames as $Value)
        {
          if ($Value!="."  AND $Value!=".." AND !is_dir($Value))
          // if (!is_dir($Value)) // this doesn't work because all dir's are symlinks!!!
          {
	    echo "<tr style=\"background-color: rgb(255,255,255);\" onmouseover=\"this.style.backgroundColor='#C4C4C4'\" onmouseout=\"this.style.backgroundColor='#FFFFFF'\" bgcolor=\"#FFFFFF\"><td>\n";
	    if ($Settings[Status] == "start")  // select the masterfile
	    {
             echo "<input type=\"radio\" name=\"MasterFile\" value=\"" .$Value ."\">\n";
	    }
	    if ($Settings[Status] == "master")  // masterfile selected, select files for aligning
	    {
             echo "<input type=\"checkbox\" name=\"Files[]\" value=\"" .$Value ."\">\n";
	    }
	    echo $Value;
	    echo "</td></tr>";
	    $cbCounter++;
          } 
        }
        echo "</table>";
      echo "</td></tr>";  
    echo "</table>";  
    switch ($Settings[Status]) {
      case "start" : $Settings[Status] = "master";
        break;
      case "master" : $Settings[Status] = "parameters";
        break;
    }
    $Settings[dir] = $DataDir;
    echo "<input type=\"hidden\" name=\"Settings\" value=\"" .base64_encode(serialize($Settings)) ."\">";

  echo "</form>";
  // HTML doesn't end here but below
}  // end if not $Settings[Status] parameters
else
{
//  print_r($Settings[Files]);
  $Settings[Status] = "startalign";
  $Settings[ShowStop] = 0; 
  echo "<b>Select paramaters for alignment :</b> ";
  echo "<form method=\"post\" action=\"msalign.php\">";
  echo "<table>";
  echo "<tr><td>";
    echo "Mass measurement error [ppm] :</td><td>";
    echo "<input type=\"text\" name=\"MassError\" value=\"15\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"The mass error in ppm's.\">Info</abbr></td></tr>";
  echo "<tr><td>";
    echo "Standard deviation retention time [scans] :<br>Set to 0 to determine this automatically</td><td>";
    echo "<input type=\"text\" name=\"TimeError\" value=\"10\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"Standard deviation of the retention time in scans.\">Info</abbr></td></tr>";
  echo "<tr><td>";
    echo "Number of features :</td><td>";
    echo "<input type=\"text\" name=\"Features\" value=\"400\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"For alignment at least 100 features are needed. The more complex the sample, the more features are nescessary.\">Info</abbr></td></tr>";
  echo "<tr><td>";
    echo "Start at scan number :</td><td>";
    echo "<input type=\"text\" name=\"StartScan\" value=\"300\">";
  echo "</td></tr>";
  echo "<tr><td>";
    echo "End at scan number :</td><td>";
    echo "<input type=\"text\" name=\"EndScan\" value=\"2100\">";
  echo "</td></tr>";
  echo "<tr><td>";
    echo "Endpoint X :</td><td>";
    echo "<input type=\"text\" name=\"EndX\" value=\"2800\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"Set the endpoint approximately 500 - 1000 higher then the end scan.\">Info</abbr></td></tr>";
  echo "<tr><td>";
    echo "Endpoint Y :</td><td>";
    echo "<input type=\"text\" name=\"EndY\" value=\"2800\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"Set the endpoint approximately 500 - 1000 higher then the end scan.\">Info</abbr></td></tr>";
  echo "<tr><td>";
    echo "Costs per breakpoint :</td><td>";
    echo "<input type=\"text\" name=\"CostsPB\" value=\"0.3\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"The costs per breakpoint determines if a breakpoint is kept or discarded. The lower the value the more breakpoints are kept. Default = 0.5\">Info</abbr></td></tr>";
  echo "<tr><td>";
    echo "Fitness :</td><td>";
    echo "<input type=\"text\" name=\"Fitness\" value=\"1\">";
  echo "</td><td><abbr style=\"border-bottom: dotted 0.05em; cursor: help;\" title=\"The fitness determines.... Default = 1\">Info</abbr></td></tr>";
  echo "<tr><td colspan=\"3\">";
    echo "<input type=\"Submit\" name=\"Submit\" value=\"Submit\"></td></tr>";
  echo "</td></tr>";
  echo "</table>";
  echo "<input type=\"hidden\" name=\"Settings\" value=\"" .base64_encode(serialize($Settings)) ."\">";
  echo "</form>";
  echo "<p><b>Selected masterfile :</b><br>" .$Settings[MasterFile] ."</p>";
  echo "<p><b>Files selected for aligment :</b><br>";
  $SelectedFiles = $Settings[Files];
  foreach ($SelectedFiles as $Value)
  {
    echo $Value ."<br>";
  }
  echo "</p>";
}


//here the HTML ends
include("footer.inc.php");
echo $footer;
?>
