# MRC/MRCS file handling in a simple and easy package - dpf-mrcfile

When working with MRC files, there are many libraries that you can work with.  Some are
faster than this one, some are slower than this one.  This library was not built out of
speed, or to be novel.  It was simply built, because I wanted a library that was easy
to interpret, easy to modify, and allowed for fast MRC file modifications.


```python

test_mrc = "test/EMD-3001.map"
m = MRC()
m.read_mrc(test_mrc)
print(m.header_as_string())
```
