#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys

def send_email(mail,PAGE_ADDRESS,jobID,date):
  
  import smtplib
  from email.mime.multipart import MIMEMultipart
  from email.mime.text import MIMEText

  # Create message container - the correct MIME type is multipart/alternative.
  msg = MIMEMultipart('alternative')
  msg['Subject'] = "MSAvis: job %s completed" % jobID
  msg['From'] = ""
  msg['To'] = mail

  # Create the body of the message (a plain-text and an HTML version).
  text = 'Your job %s is complete (from %s).<p>You can fing your image <a href=%s/msa_vis/result/%s>here</a> or download below:' % (jobID, date, PAGE_ADDRESS, jobID)
  html = """\
  <html>
    <head>Your job %s is completed (from %s).</head>
    <body>
      <p>You can fing your image <a href=%smsa_vis/result/%s>here</a>.
    </body>
  </html>
  """ % (jobID, date, PAGE_ADDRESS, jobID)

  # Record the MIME types of both parts - text/plain and text/html.
  part1 = MIMEText(text, 'plain')
  part2 = MIMEText(html, 'html')

  # Attach parts into message container.
  # According to RFC 2046, the last part of a multipart message, in this case
  # the HTML message, is best and preferred.
  msg.attach(part1)
  msg.attach(part2)

  # Send the message via SMTP server.
  port = 465
  host = 'smtp.gmail.com'
  user = ""
  password = ""
  s=smtplib.SMTP_SSL(host,port)
  s.login(user, password)

  s.sendmail(msg['From'],msg['To'],msg.as_string())
  s.quit()
 


send_email(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]+" "+sys.argv[5])  
  