from django.db import models
# import django.db.models.JSONField
from django.db.models import JSONField


# Create your models here.
class Item(models.Model):
    species = models.CharField(max_length=200)
    gene_symbol = models.CharField("Gene name", max_length=200)
    gene_id = models.CharField(max_length=200)
    start = models.CharField("Start of interval")
    end = models.CharField("End of interval")
    snp_id = models.CharField("SNP ID")
    snp_pos = models.CharField("SNP localization")
    data = JSONField()
    article_title = models.CharField(max_length=200)
    article_url = models.CharField(max_length=200)

    # pub_date = models.DateTimeField("date published")
    # votes = models.IntegerField(default=0)

    # start = models.IntegerField("Start of interval", default=0)
    # end = models.IntegerField("End of interval", default=0)

    # start = models.FloatField("Start of interval")
    # end = models.FloatField("End of interval")


# species
# gene_name
# gene_id
# start
# end
# snp_id
# snp_pos
# data
# article_title
# article_url



# Gene name
# Species
# Gene ID
# Start of interval
# End of interval
# SNP ID
# SNP localization
# Positive selection metrics

