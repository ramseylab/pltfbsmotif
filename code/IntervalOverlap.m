function f=IntervalOverlap(s1, e1, s2, e2)

if e1 < s1
  error 'invalid data';
end

if e2 < s2
  error 'invalid data';
end

if s1 < s2
  if e1 < s2
    f = 0;
  else
    f = (min(e2,e1)-s2)/(e1-s1);
  end
else
  if e2 < s1
    f = 0;
  else
    f = (min(e2,e1)-s1)/(e1-s1);
  end
end
