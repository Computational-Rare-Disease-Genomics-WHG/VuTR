# nginx/Dockerfile 

# configuration file from https://github.com/docker/awesome-compose/blob/master/nginx-wsgi-flask/nginx/Dockerfile
FROM nginx:1.19.7-alpine

# Add nginx.conf to container
COPY --chown=nginx:nginx nginx.conf /etc/nginx/nginx.conf
COPY --chown=nginx:nginx start.sh /app/start.sh

# set workdir
WORKDIR /app

# permissions and nginx user for tightened security
RUN chown -R nginx:nginx /app && chmod -R 755 /app && \
        chown -R nginx:nginx /var/cache/nginx && \
        chown -R nginx:nginx /var/log/nginx && \
        chmod -R 755 /var/log/nginx; \
        chown -R nginx:nginx /etc/nginx/conf.d
RUN touch /var/run/nginx.pid && chown -R nginx:nginx /var/run/nginx.pid

USER nginx

CMD ["nginx", "-g", "'daemon off;'"]
