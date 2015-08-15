#! /usr/bin/python

import cgitb
import cgi
import shlex
import sys, os
sys.stderr = sys.stdout

def header():
  print "Content-Type: text/html"
  print
  print "<html>"
  print "<head><title>sginfo gateway</title></head>"
  print "<body><h2>SgInfo - Space group Info</h2><hr>"

def please_visit_cctbx():
  print '<strong><a href="http://cci.lbl.gov/cctbx/"'
  print ">Please visit our more recent cctbx services</a></strong>"

def tail():
  print "<hr>"
  please_visit_cctbx()
  print "<hr>"
  print "<address>"
  print "Ralf W. Grosse-Kunstleve"
  print "&lt;<a href=\"mailto:rwgkio+sginfo@gmail.com\""
  print "                    >rwgkio+sginfo@gmail.com</a>&gt;"
  print "</address>"
  print "</body>"
  print "</html>"

def show_string(s):
  if (s is None): return None
  if (s.find('"') < 0): return '"'+s+'"'
  if (s.find("'") < 0): return "'"+s+"'"
  return '"'+s.replace('"','\\"')+'"'

def quote_arg(arg):
  try: flds = shlex.split(arg)
  except: flds = []
  if (len(flds) == 1): return arg
  return show_string(arg)

def str_to_html(s):
  return (s
    .replace('&', '&amp;')
    .replace('<', '&lt;')
    .replace('>', '&gt;')
    .replace('"', '&quot;'))

def input_form(args):
  print '<form method="post" action="sginfo-query.cgi">'
  print '<code><strong>sginfo </strong></code><input'
  print 'name="argl" size=60'
  print 'value="%s"><p>' % str_to_html(
    " ".join([quote_arg(arg) for arg in args]))
  print '<input type="submit" value="Run sginfo">'
  print '<input type="reset"  value="Reset">'
  print '<strong><a href="/sginfo/">'
  print 'On-line Documentation</a></strong>'
  print '</form>'
  print '<hr>'
  please_visit_cctbx()
  print '<hr>'

def run():
  cgitb.enable()
  header()
  form = cgi.FieldStorage()
  try: argl = form['argl'].value
  except KeyError: args = []
  else: args = shlex.split(argl)
  input_form(args)
  args.insert(0, "sginfo")
  os.environ["PATH"] = os.getcwd()
  child_stdin, child_stdout_and_stderr = os.popen4(args)
  child_stdin.close()
  output = str_to_html(child_stdout_and_stderr.read())
  print "<pre>"
  sys.stdout.write(output)
  print "</pre>"
  tail()

if (__name__ == "__main__"):
  run()
