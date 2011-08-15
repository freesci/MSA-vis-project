from django.conf.urls.defaults import *

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Example:
     #(r'^msa_vis/', include('msa_vis.foo.urls'))
     (r'^msa_vis/$', 'msa_vis.msa_vis_app.views.first_page'),
     (r'^msa_vis/edit/$', 'msa_vis.msa_vis_app.views.second_page'),
     (r'^chart/simple.png$', 'msa_vis_app.chart.simple')

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # (r'^admin/', include(admin.site.urls)),
)
