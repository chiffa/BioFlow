import logging
from smtplib import SMTP
from datetime import datetime
from email.message import EmailMessage

from bioflow.configs.main_configs import smtp_logging, smtp_logging_parameters

# define a formatter
formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')


mime_message = EmailMessage()
mime_message['From'] = smtp_logging_parameters['local_mail_account']
mime_message['To'] = smtp_logging_parameters['reporting_target_mail']


def get_smtp_logger(logger_name):
    """
    Returns a properly configured logger object
    :param logger_name: name of the logger object
    """
    _logger = logging.getLogger(logger_name)
    _logger.setLevel(logging.INFO)

    if smtp_logging:
        mail_handler = logging.handlers.SMTPHandler(
            mailhost=smtp_logging_parameters['local_host'],
            fromaddr=smtp_logging_parameters['local_mail_account'],
            toaddrs=smtp_logging_parameters['reporting_target_mail'],
            subject="BioFlow runtime error"
        )

        mail_handler.setLevel(logging.ERROR)
        _logger.addHandler(mail_handler)

    return _logger


def successfully_completed(start_date, start_device):

    with SMTP(host=smtp_logging_parameters['local_host']) as smtp_server:
        mime_message['Subject'] = 'Run started on %s and %s has completed' % (start_date.isoformat(), start_device)
        mime_message.set_content('Run has completed successfully after %s minutes' % ((datetime.now() -
                   start_date).total_seconds()/60.))

        smtp_server.send_message(mime_message)


def smtp_error_bail_out():

    with SMTP(host=smtp_logging_parameters['local_host']) as smtp_server:
        mime_message['Subject'] = 'SMTPHandler error bail-out'
        mime_message.set_content("There was an error in the code and the logger's SMPTHandler bailed out.")
        print(smtp_server.noop())
        smtp_server.send_message(mime_message)