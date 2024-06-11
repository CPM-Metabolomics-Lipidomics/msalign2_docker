<?PHP
include("functions.inc.php");
$beginheader = "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\">\n";
$beginheader .= "<HTML>\n";
$beginheader .= "<HEAD>\n";
$beginheader .= "<META HTTP-EQUIV=\"content-type\" CONTENT=\"text/html; charset=UTF-8\">\n";
$beginheader .= "<TITLE>LUMC - PARA - MS Alignment</TITLE>\n";
$scriptheader = "<style type=\"text/css\">\n";
$scriptheader .= "<!--";
$scriptheader .= "@import url(\"./styles/style.css\");\n";
$scriptheader .= "-->";
$scriptheader .= "</style>\n";
$scriptheader .= "<script language=\"javascript\">\n<!--\n";
$scriptheader .= "function SetAllCheckBoxes(FormName, FieldName, CheckValue)\n";
$scriptheader .= "{\n";
$scriptheader .= "  if(!document.forms[FormName])\n  return;\n";
$scriptheader .= "  var objCheckBoxes = document.forms[FormName].elements[FieldName];\n";
$scriptheader .= "  if(!objCheckBoxes)\n  return;\n";
$scriptheader .= "  var countCheckBoxes = objCheckBoxes.length;\n";
$scriptheader .= "  if(!countCheckBoxes)\n  objCheckBoxes.checked = CheckValue;\n";
$scriptheader .= "  else\n";
$scriptheader .= "    for (var i=0;i<countCheckBoxes;i++)\n  objCheckBoxes[i].checked = CheckValue;\n";
$scriptheader .= "}\n";
$scriptheader .= "-->\n</script>\n";

//-------------------------- get the server load ----------------------------------------------------
  $GetLoad = "total_load.sh";
  exec ($GetLoad, $ServerLoad);
  $iServerLoad = 1-(str_replace("%","",$ServerLoad[0])/100);     // 0 is read, 1 is green so correct by substracting the server load from 1
//------------------------- end server load ---------------------------------------------------------

$bodyheader = "</HEAD>\n";
$bodyheader .= "<BODY class=\"main\">\n";
$bodyheader .= "<table width=\"100%\" class=\"main\">\n";
$bodyheader .= "<tr><td class=\"menu\"><a href=\"index.php\"><img src=\"./images/lumc_logo.jpg\" align=\"top\" style=\"border-style: none;\"></a>\n";
$bodyheader .= "<a class=\"helplink\" href=\"help.php\">Help</a>";
$bodyheader .= "</td><td align=\"center\" width=\"150\" style=\"background-color: rgb(" .CalcColor($iServerLoad) .");\">"; 
$bodyheader .= "<b>Server load : " .sprintf("%0.1f", $ServerLoad[0]) ."&#37;</b>";
$bodyheader .= "</td></tr>\n";
$bodyheader .= "</table>\n";

$nobodyheader = "</head><body class=\"picture\">";
?>
