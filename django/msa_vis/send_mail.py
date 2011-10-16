import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import settings


def send_email_to(mail,settingsPAGE_ADRESS,jobID,date):

    # Create message container - the correct MIME type is multipart/alternative.
    msg = MIMEMultipart('alternative')
    msg['Subject'] = "MSAvis"
    msg['From'] = "noreply"
    msg['To'] = mail

    # Create the body of the message (a plain-text and an HTML version).
    text = 'Your job %s is complete (from %s).<p>You can fing your image <a href=%s/msa_vis/result/%s>here</a>' % (jobID, date, settingsPAGE_ADRESS, jobID)
    html = """\
    <html>
      <head>Your job %s is complete (from %s).</head>
      <body>
	<p>You can fing your image <a href=%smsa_vis/result/%s>here</a> </p>
      </body>
    </html>
    """ % (jobID, date, settingsPAGE_ADRESS, jobID)

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