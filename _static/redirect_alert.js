// Modified minimalist version of https://github.com/usnistgov/leaveNotice/blob/nist-pages/js/jquery.leaveNotice-nist.min.js

(function($){$.fn.leaveNotice=function(opt){
	var defaults={};
	var options=$.extend(defaults,opt);

	return this.each(function(){el=$(this);
		// get the url
		var url=el.attr("href");
		
		// create a shortened version of the url
		var ulen=options.displayUrlLength;
		if(url.length>=ulen){var suffix="..."}else{var suffix=""}
		var shortUrl=url.substr(0,ulen)+suffix;
		
		// display alert if clicked
		el.click(function(){
		        if (confirm("Thank you for visiting NIST.\n\nWe have provided a link to this site because it has information that may be of interest to our users. NIST does not necessarily endorse the views expressed or the facts presented on this site. Further, NIST does not endorse any commercial products that may be advertised or available on this site.\n\nClick OK to be directed to: "+shortUrl)) {window.location = url};
			return false})
		})
	};
})(jQuery);
