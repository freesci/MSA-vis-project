from django.conf.urls.defaults import *

# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('',
    # Example:
     #(r'^msa_vis/', include('msa_vis.foo.urls'))
     (r'^msa_vis/$', 'msa_vis.msa_vis_app.views.first_page'),
     (r'^msa_vis/edit/$', 'msa_vis.msa_vis_app.views.second_page'),
     (r'^chart/simple.png$', 'msa_vis_app.chart.simple'),
     (r'^msa_vis/error1/$', 'msa_vis.msa_vis_app.views.third_page'),
     (r'^chart_download/simple.png$', 'msa_vis_app.chart_download.simple')

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # (r'^admin/', include(admin.site.urls)),
)

from django.conf import settings
if settings.DEBUG: 
    urlpatterns += patterns('',
        (r'^static/(?P<path>.*)$', 'django.views.static.serve', {'document_root': settings.MEDIA_ROOT}),
    )
