<?PHP
if ((isset($_POST['Settings'])) || (isset($_GET['Settings'])))
{
  if (isset($_POST['Settings']))
  {
    $Settings = unserialize(base64_decode($_POST['Settings']));
  }
  else
  {
    $Settings = unserialize(base64_decode($_GET['Settings']));
  }
  $filenames = $Settings[Files];
  $NumberFiles = count($filenames);
  $allClear = $NumberFiles + 2;  // +2 is for temp file and running file
  foreach ($filenames as $value)
  { 
    $nametmp = str_ireplace(".mzxml", "_out.txt", $value);
    $do = unlink($nametmp);
    if ($do == 1)
    {
      $allClear--;
    }
    else
    {
      echo "There was an error!";
    }
  }
  $do = unlink($Settings[tmpFile]);
  if ($do == 1)
  {
    $allClear--;
  }
  else
  {
    echo "There was an error!";
  }
  $do = unlink($Settings[runningFile]);
  if ($do == 1)
  {
    $allClear--;
  }
  else
  {
    echo "There was an error!";
  }

  if ($allClear==0)
  {
    $NewLocation = "location: index.php";
    header($NewLocation);
  }
  else
  {
    echo "Not all files where deleted!!";
  }


}
?>
