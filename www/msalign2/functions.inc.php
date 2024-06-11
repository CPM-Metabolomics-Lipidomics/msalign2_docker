<?php

function utime() {
  $time = explode(" ", microtime());
  $usec = (double)$time[0];
  $sec = (double)$time[1];
  return $sec + $usec;
}

// create hex color from rgb code
function rgb2hex($red,$green,$blue) {
  $red = dechex($red);
  $green = dechex($green);
  $blue = dechex($blue);
  return "#" .strtoupper($red) .strtoupper($green) .strtoupper($blue);
}
// calculate color according to Quality Factor
function CalcColor($Quality) {
  // color should go from red to yellow to green -> bad to medium to good
  // red = 255,0,0 / yellow = 255,255,0 / green = 0,255,0
  if ($Quality < 0.5)
  {
    $red = 255;
    $green = round((2*$Quality) * 255);
  }
  else
  {
    $red = round(2*(1-$Quality) * 255);
    $green = 255;
  }
  $blue = 0;
//  return rgb2hex($red, $green, $blue);
  $color = $red ."," .$green ."," .$blue;
  return $color;
}

function FindString ($String, $Findme){
  $pos = strpos($String, $Findme);
  return $pos;
}

function OddEven ($x){
  if ($x & 1) return "odd";
    else return "even";
}


?>
