<?php
include("header.inc.php");
echo $beginheader;
echo $scriptheader;
echo $bodyheader;

if (isset($_GET['Settings']))
{
  $Settings= unserialize(base64_decode($_GET['Settings']));

  echo "<a href=\"index.php\">HOME</a>";
  echo "<hr>";

if ($Settings[Status] == "align")
{
  $StatusPage = "";

  $FilesOut = $Settings[Files];
  $StatusFiles = "<table border=\"1\">";
  $NumberFiles = count($FilesOut);
  foreach ($FilesOut as $value)
  {
    $StatusFiles .= "<tr><td><b>" .$value ."</b><br>";
    $FileNameTXT = str_ireplace(".mzxml", "_out.txt", $value);
  
    $checkfile = FALSE;
    $i = 0;
    while(($checkfile == FALSE) && ($i <= 5))
    {
      usleep(50000);
      if (file_exists($FileNameTXT))
      {
        $info = file($FileNameTXT);
        foreach ($info as $line)
        {
          $StatusFiles .= $line ."<br>";
	  if (FindString($line, "Finished") !== FALSE)
	  {
            $NumberFiles--;
	    $StatusFiles .= "<a href=\"" .$value .".alignment.png\">Check alignment</a>&nbsp;&nbsp;&nbsp;<a href=\"" .$value .".alignment\">Check mass pairs</a>";
	  }
	  if (FindString($line, "aborted") !== FALSE)
	  {
            $NumberFiles--;
	  }
        } // end foreach $info
        $checkfile = TRUE;
      }  // end if is_readable
      else
      {
	$i = 6;
      }
      $i++;
    } // end of while $checkfile loop
    $StatusFiles .= "</td></tr>";
  }  // end foreach $FilesOut

  $StatusFiles .= "</table>";
  if ($NumberFiles == 0)     // show the continue button if all files are done
  {
    $Settings[ShowStop] = 1;
  }
  if ($Settings[ShowStop] == 0)
  {
    $StatusPage = "<a href=\"status.php?Settings=" .base64_encode(serialize($Settings)) ."\">Refresh page</a><br>";
  }
  if ($Settings[ShowStop] == 1)     // show the continue button if all files are done
  {
    $StatusPage .= "<form action=\"cleartmp.php\" method=\"post\">";
    $StatusPage .= "<input type=\"Submit\" name=\"Submit\" value=\"Remove temporary files\">";
    $StatusPage .= "<input type=\"hidden\" name=\"Settings\" value=\"" .base64_encode(serialize($Settings)) ."\">";
    $StatusPage .= "</form>";
  }
  echo $StatusPage .$StatusFiles;
}  // end status==align
if ($Settings[Status]=="startalign")
{
  $Settings[Status] = "align";
  echo "<a href=\"status.php?Settings=" .base64_encode(serialize($Settings)) ."\">Refresh page</a>";
  echo "<p>Calculations running. Refresh this page to see the progress of your calculations.</p>";
}
}
?>
