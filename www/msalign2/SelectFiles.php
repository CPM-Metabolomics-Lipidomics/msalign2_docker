<?PHP

if (!isset($_POST['Submit']))
{
  $DataDir = $_GET['dir'];
  if ($DataDir=="") $DataDir = "data";

  echo $DataDir ."<br><br>";

  $ListDirs = "ls -H ./data/";
  exec ($ListDirs, $Dirs);

  echo "<a href=\"SelectFiles.php?dir=Data\">Data</a>";
  foreach ($Dirs as $Value)
  {
    echo "&nbsp;&nbsp;<b>/</b>&nbsp;&nbsp;<a href=\"SelectFiles.php?dir=data/".$Value ."\">" .$Value ."</a>\n";
  }
  echo "<BR><BR>\n";

  // formulate the command, make sure to add > for redirecting output to a file
  $ListFiles = "ls -al ./" .$DataDir ."/";
  
  
  // executes the command and writes the input to a file
  exec($ListFiles, $FileNames);

  echo "Please select files to submit";
  echo "<form method=\"post\" action=\"test.php\">";
  foreach ($FileNames as $Value)
  {
    if ($Value!="."  AND $Value!=".." AND !is_dir($Value))
    // if (!is_dir($Value)) // this doesn't work because all dir's are symlinks!!!
    {
      echo "<input type=\"checkbox\" name=\"file[ ]\" value=\"" .$Value ."\">\n";
      echo $Value ."<br>\n";
    } 
  }
  echo "<br><input type=\"Submit\" name=\"Submit\" value=\"Submit\">";
  echo "</form>";
}
else
{
  echo "Selected files are :<br>";
  $filenames = $_POST['file'];
  foreach ($filenames as $value)
  {
    echo $value ."<br>";
 }
}

?>
