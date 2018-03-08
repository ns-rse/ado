program define _gamenu_error
  version 7.0
  global GA_menu_error "`*'"
  window control static GA_menu_error 10 5 180 10 center 
  global GA_menu_OK "exit 3001"
  window control button "OK"  80 20 40 20 GA_menu_OK default
  cap noisily window dialog "Error" 100 100 200 45
  macro drop GA_menu_*
  end
