let m4p = (function(d3){

  // Calculate the percentiles of a numeric vector
  let percentile = function(d,p){
		d.sort(d3.ascending);
		// console.log(d);
		index=Math.round((d.length*p)/100)-1;
		// console.log(index);
		// console.log(d[index]);
		return(d[index]);
	};

  return {
    percentile: percentile
  };
})(d3);
