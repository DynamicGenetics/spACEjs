// m4p module for spACE
// Copyright (C) Oliver Davis 2020

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version. See https://www.gnu.org/licenses/

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
