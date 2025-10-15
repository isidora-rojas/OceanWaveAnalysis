function [year,month,day,hours,minutes,seconds,msecs]=str2dat(Str,Format);
% str2dat    : convert string to date
%
% DESCRIPTION: [year,month,day,hours,minutes,seconds,msecs]=str2dat(Str,Format);
%
% PARAMETERS :
%    year,month,day,hours, minutes and seconds are trivial
%    msecs indicates milliseconds
%    str       the string to convert (see examples)
%    format    any combination of characters, of which the
%              following character series are used for formatting:
%              yy   mod(year,100)
%              yyyy the 4 digit year
%              mm   month number
%              mmm  3 character month name
%              dd   day number
%              hh   hours
%              tt   minutes
%              ss   seconds
%              msc  milliseconds
%
% EXAMPLES   :
% [1990,10,12,0,0,0,0] = str2dat('12-10-1990','dd-mmm-yyyy')
% [1990,10,12,0,0,0,0] = str2dat('901012,'yymmdd')
% [0,0,0,17,21,0,0]    = str2dat('17:21h','hh:tth');
% [0,0,0,2,5,23,100]   = str2dat('020523_100','hhttss_msc');
%

ValidForm  = ['yy';'mm';'dd';'hh';'tt';'ss'];

MonthNames = ['JAN';'FEB';'MAR';'APR';'MAY';'JUN';...
              'JUL';'AUG';'SEP';'OCT';'NOV';'DEC' ];

ScanFormat = ['%4d';'%3s';'%2d';'%2d';'%2d';'%2d';'%2d';'%2d';'%4d'];
dat        = -ones(size(ScanFormat,1));

str        = Str;
format     = Format;
LenForm    = length(format);

while(LenForm>0),

   if ( strcmp(format(1:min(4,LenForm)),'yyyy') ),
      % format like 1990
      n = 1;
   elseif ( strcmp(format(1:min(3,LenForm)),'mmm') ),
      % format like OCT
      n = 2;
   elseif ( strcmp(format(2:min(4,LenForm)),'msc') ),
      % format like 100
      n = 9;
   else
      i = strmatch(format(1:min(2,LenForm)),ValidForm,'exact');
      if length(i)==1
         n = i + 2;
      else
         n = 10;
      end
   end

%  Month name ?
   if n==2
      [MonthStr,no,err,next]=sscanf(str,ScanFormat(n,:),1);
      i = strmatch(upper(MonthStr),MonthNames);
      if length(i)~=1
         error('invalid month');
         i = -1;
      end
      dat(n) = i;

%  Character not in ValidForm ?
   elseif n==10
      next = 2;

%  Scan number
   else
      if str(1)==' '; str(1) = '0'; end;
      if str(1)=='_'; str(1) = '0'; end;
      [dat(n),no,err,next]=sscanf(str,ScanFormat(n,:),1);
   end

   str    = str(next:length(str));
   format = format(next:LenForm);
   LenForm= LenForm - next + 1;

end; %while

if dat(3)>=0; year = dat(3)+1900; else year = -1; end;
month   = dat(4);
day     = dat(5);
hours   = dat(6);
minutes = dat(7);
seconds = dat(8);
msecs   = dat(9);
if dat(1)>=0; year  = dat(1); end;
if dat(2)>=0; month = dat(2); end;
