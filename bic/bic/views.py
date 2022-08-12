from django.http import HttpResponse, JsonResponse
from django.template import loader
from django.shortcuts import render
from analyses.models import Statistics

def index(request):
    title = "BIC"
    context = {"title": title}
    return render(request, 'index.html', context=context)

def statistics(request):
    title = "Statistics"

    statistics_query_set = Statistics.objects.order_by('cancer_type')
    context = {"title": title, "data": statistics_query_set}

    return render(request, 'statistics.html', context=context)
    


def documents(request):
    title = "Documents"
    context = {"title": title}
    return render(request, 'documents.html', context=context)

def downloads(request):
    title = "Downloads"
    context = {"title": title}
    return render(request, 'downloads.html', context=context)

# def contact(request):
#     title = "Contact"
#     context = {"title": title}
#     return render(request, 'contact.html', context=context)


