$(document).ready(function(){
  // Mark external (non-nist.gov) A tags with class "external"
  //If the adress start with https and ends with nist.gov
  var re_nist = new RegExp('^https?:\/\/((^\/)*\.)*nist\\.gov(\/|$)');
  //Regex to find address that start with https
  var re_absolute_address = new RegExp('^((https?:)?\/\/)');
  $("a").each(function(){
    var url=$(this).attr('href');
    if(re_nist.test(url) || !re_absolute_address.test(url)){
      $(this).addClass('local');
    }else{
      $(this).addClass('external');
    }
  });
  // Add leaveNotice to external A elements
  $('a.external').leaveNotice({
    siteName: 'pages.nist.gov',
  });
});
