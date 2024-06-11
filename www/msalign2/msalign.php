<?PHP

if (!isset($_POST['Settings']))
{
  echo "Nothing processed!!";
}
else
{
// go to statuspage
  $Settings = unserialize(base64_decode($_POST['Settings'])); 
// write the filenames to a file 1 on each line
  $UniqID = uniqid();
  $NameTempFile = "./temp/" .$UniqID .".txt";
  $NameRunningFile = "./temp/running" .$UniqID .".txt";
  $Settings[tmpFile] = $NameTempFile;
  $Settings[runningFile] = $NameRunningFile;
  $Settings[UniqID] = $UniqID;

  $NewHeader = "Location: status.php?Settings=" .base64_encode(serialize($Settings));
  header($NewHeader);

  $filenames = $Settings[Files];
  $features = $_POST['Features'];
  $CostsPB = $_POST['CostsPB'];
  $fitness = $_POST['Fitness'];
  $mass_error = $_POST['MassError'];
  $startscan = $_POST['StartScan'];
  $endscan = $_POST['EndScan'];
  $timewindow = $_POST['TimeError'];
  $masterfile = $Settings[MasterFile];
  $EndX = $_POST['EndX'];
  $EndY = $_POST['EndY'];

// create file with all the set parameters
  $ParamFileName = $Settings[dir] ."/AlignmentParameters.txt";
  $ParamFile = fopen($ParamFileName, 'w');
  $ParamText = "Mass measurement error : " .$mass_error ."\n";
  fwrite($ParamFile, $ParamText);
  $ParamText = "Standard deviation : " .$timewindow ."\n";
  fwrite($ParamFile, $ParamText);
  $ParamText = "Number of features : " .$features ."\n";
  fwrite($ParamFile, $ParamText);
  $ParamText = "Scan range : " .$startscan ." - " .$endscan ."\n";
  fwrite($ParamFile, $ParamText);
  $ParamText = "Endpoint X, Y : " .$EndX .", " .$EndY ."\n";
  fwrite($ParamFile, $ParamText);
  $ParamText = "Costs per breakpoints : " .$CostsPB ."\n";
  fwrite($ParamFile, $ParamText);
  $ParamText = "Fitness : " .$fitness ."\n";
  fwrite($ParamFile, $ParamText);
  fclose($ParamFile);


  $MyTempFile = fopen($NameTempFile, 'w');
  foreach ($filenames as $Value)
  {
    $stringData = $Value ."\n";
    fwrite($MyTempFile, $stringData);
  }
  fclose($MyTempFile);   // close the temporary file
  $command = "./MSAlign.sh " .$NameTempFile ." " .$masterfile ." " .$mass_error ." " .$timewindow ." " .$features ." " .$CostsPB ." " .$fitness ." " .$EndX ." " .$EndY ." " .$startscan ." " .$endscan ." > " .$NameRunningFile ." &";
  exec($command);
}

?>
